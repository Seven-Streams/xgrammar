import itertools
import json
from typing import Optional, Sequence

import pytest

import xgrammar as xgr
from xgrammar.testing import _get_masked_tokens_from_bitmask


def _compile_lark(
    grammar: str, tokenizer_info: Optional[xgr.TokenizerInfo] = None
) -> xgr.CompiledGrammar:
    tokenizer_info = tokenizer_info or xgr.TokenizerInfo([])
    grammar_obj = xgr.Grammar.from_lark(grammar, tokenizer_info=tokenizer_info)
    return xgr.GrammarCompiler(tokenizer_info, cache_enabled=False).compile_grammar(grammar_obj)


def _matches_string(compiled: xgr.CompiledGrammar, value: str) -> bool:
    matcher = xgr.GrammarMatcher(compiled, terminate_without_stop_token=True)
    return matcher.accept_string(value) and matcher.is_terminated()


def _assert_grammar_language(
    grammar: xgr.Grammar,
    accepted: Sequence[str],
    rejected: Sequence[str],
    tokenizer_info: Optional[xgr.TokenizerInfo] = None,
) -> None:
    tokenizer_info = tokenizer_info or xgr.TokenizerInfo([])
    compiled = xgr.GrammarCompiler(tokenizer_info, cache_enabled=False).compile_grammar(grammar)
    for value in accepted:
        assert _matches_string(compiled, value), value
    for value in rejected:
        assert not _matches_string(compiled, value), value


def _assert_language(
    grammar: str,
    accepted: Sequence[str],
    rejected: Sequence[str],
    tokenizer_info: Optional[xgr.TokenizerInfo] = None,
) -> None:
    grammar_obj = xgr.Grammar.from_lark(grammar, tokenizer_info=tokenizer_info)
    _assert_grammar_language(grammar_obj, accepted, rejected, tokenizer_info)


def _matches_token_sequence(compiled: xgr.CompiledGrammar, token_ids: Sequence[int]) -> bool:
    matcher = xgr.GrammarMatcher(compiled, terminate_without_stop_token=True)
    for token_id in token_ids:
        if not matcher.accept_token(token_id):
            return False
    return matcher.is_terminated()


def _assert_token_language(
    grammar: str,
    tokenizer_info: xgr.TokenizerInfo,
    accepted: Sequence[Sequence[int]],
    rejected: Sequence[Sequence[int]],
) -> None:
    compiled = _compile_lark(grammar, tokenizer_info)
    for token_ids in accepted:
        assert _matches_token_sequence(compiled, token_ids), token_ids
    for token_ids in rejected:
        assert not _matches_token_sequence(compiled, token_ids), token_ids


def _assert_lark_error(
    grammar: str, message: str, tokenizer_info: Optional[xgr.TokenizerInfo] = None
) -> str:
    with pytest.raises(RuntimeError) as exc_info:
        xgr.Grammar.from_lark(grammar, tokenizer_info=tokenizer_info)
    error = str(exc_info.value)
    assert message in error
    assert "Lark error at line " in error
    assert ", column " in error
    return error


@pytest.mark.parametrize(
    "grammar, accepted, rejected",
    [
        pytest.param(
            'start: "a" "b" | "c" ("d" | "e")',
            ["ab", "cd", "ce"],
            ["", "a", "c", "ade"],
            id="sequence-choice-precedence",
        ),
        pytest.param(
            'start: ("a" | "b") ["c"]',
            ["a", "b", "ac", "bc"],
            ["", "c", "abc"],
            id="groups-and-optional-group",
        ),
        pytest.param(
            'start: | "a" | "b"', ["", "a", "b"], ["ab", "c"], id="empty-left-alternative"
        ),
        pytest.param('start: "a" |', ["", "a"], ["aa", "b"], id="empty-right-alternative"),
        pytest.param("start:", [""], ["a", " "], id="empty-rule"),
        pytest.param('start: "" "a" ""', ["a"], ["", "aa"], id="empty-literals"),
        pytest.param(
            '?item: "a"\n!suffix: "b"\nstart: item suffix',
            ["ab"],
            ["", "a", "b"],
            id="lark-rule-prefixes",
        ),
        pytest.param(
            'start: _item _TOKEN\n_item: "a"\n_TOKEN: "b"',
            ["ab"],
            ["a", "b", "_item_TOKEN"],
            id="hidden-rule-and-terminal-names",
        ),
        pytest.param(
            'start: "a" -> first\n     | "b" -> second',
            ["a", "b"],
            ["", "first", "second"],
            id="alternative-aliases",
        ),
        pytest.param(
            'start: value\nvalue: "x" | "(" value ")" | "[" values? "]"\nvalues: value ("," value)*',
            ["x", "(x)", "((x))", "[]", "[x]", "[x,(x),[x]]"],
            ["", "()", "[", "[x,]", "[(])"],
            id="forward-references-and-recursion",
        ),
        pytest.param(
            'start: sequence\nsequence: "x" | "a" sequence "b"',
            ["x", "axb", "aaxbb"],
            ["", "ab", "aaxb"],
            id="recursive-sequence",
        ),
        pytest.param(
            'start: SIGNED\nSIGNED: SIGN? DIGIT+\nSIGN: "+" | "-"\nDIGIT: "0".."9"',
            ["0", "123", "+7", "-42"],
            ["", "+", "--1", "1a"],
            id="terminal-composition",
        ),
        pytest.param(
            "start: \"'foo'\" /a+/ | STRING /b+/\nSTRING: /'[^']*'/",
            ["'foo'a", "'foo'aaa", "'bar'b", "'bar'bbb", "'foo'bb"],
            ["'bar'a", "'bar'c", "foo"],
            id="literal-terminal-ambiguity",
        ),
        pytest.param(
            'start: /.../ "abc" /.../',
            ["abcabcabc", "aaaabcccc", "🔵🟠✅abc❌🟠🔵"],
            ["aaabcccc", "aaaaabcccc", "🔵🟠abc🟠🔵"],
            id="regex-dot-counts-unicode-codepoints",
        ),
        pytest.param(
            r"start: /a\/b/ /[0-9]{2,4}/",
            ["a/b12", "a/b1234"],
            ["a/b1", "a/b12345", "a-b12"],
            id="regex-escaped-delimiter-and-repeat",
        ),
        pytest.param(
            'start: "a".."z"+ "0".."9"?',
            ["a", "xyz", "hello7"],
            ["", "A", "7", "abc78"],
            id="ascii-character-ranges",
        ),
        pytest.param(
            'start: "α".."γ"+', ["α", "βγ", "γα"], ["", "δ", "a"], id="unicode-character-ranges"
        ),
        pytest.param(
            'start: "😀" "é" "中文"',
            ["😀é中文"],
            ["", "😀é", "😀e中文"],
            id="unicode-string-literals",
        ),
        pytest.param(
            r'''start: "\n" "\t" "\\" "\"" "\u03bb"''',
            ['\n\t\\"λ'],
            [r"\n\t\"λ", "\n\tλ"],
            id="json-style-string-escapes",
        ),
        pytest.param(
            r'''start: "\b" "\f" "\r"''',
            ["\b\f\r"],
            ["bfr", "\b\f\n"],
            id="control-character-string-escapes",
        ),
        pytest.param(
            'start: foo-bar FOO-BAR\nfoo-bar: "a"\nFOO-BAR: "b"',
            ["ab"],
            ["", "a-b", "foo-barFOO-BAR"],
            id="hyphenated-identifiers",
        ),
        pytest.param(
            'start: item // top-level comment\nitem: "ok" # rule comment',
            ["ok"],
            ["", "item", "ok#"],
            id="comment-styles",
        ),
    ],
)
def test_lark_core_languages(
    grammar: str, accepted: Sequence[str], rejected: Sequence[str]
) -> None:
    _assert_language(grammar, accepted, rejected)


@pytest.mark.parametrize(
    "grammar, accepted, rejected",
    [
        pytest.param('start: "a"?', ["", "a"], ["aa"], id="question"),
        pytest.param('start: "a"*', ["", "a", "aaaa"], ["b", "aaab"], id="star"),
        pytest.param('start: "a"+', ["a", "aaaa"], ["", "b"], id="plus"),
        pytest.param('start: "a"~2', ["aa"], ["", "a", "aaa"], id="tilde-exact"),
        pytest.param(
            'start: "a"~2..4', ["aa", "aaa", "aaaa"], ["", "a", "aaaaa"], id="tilde-range"
        ),
        pytest.param('start: "a"{2}', ["aa"], ["", "a", "aaa"], id="brace-exact"),
        pytest.param(
            'start: "a"{2,4}', ["aa", "aaa", "aaaa"], ["", "a", "aaaaa"], id="brace-range"
        ),
        pytest.param('start: "a"{2,}', ["aa", "aaaaaa"], ["", "a"], id="brace-open-end"),
        pytest.param('start: "a"{,2}', ["", "a", "aa"], ["aaa"], id="brace-open-start"),
        pytest.param('start: "a"{0}', [""], ["a"], id="zero-exact"),
        pytest.param('start: "a"{0,0}', [""], ["a"], id="zero-range"),
        pytest.param(
            'start: ("a" | "bc"){2,3}',
            ["aa", "abc", "bca", "bcbcbc"],
            ["", "a", "bcbc bcbc", "aaaa"],
            id="group-repeat",
        ),
        pytest.param(
            'start: ITEM{2,3}\nITEM: "x" | "y"',
            ["xx", "xy", "yyy"],
            ["", "x", "xxxx"],
            id="terminal-repeat",
        ),
        pytest.param(
            'start: item{2,3}\nitem: "x" | "(" item ")"',
            ["xx", "x(x)", "(x)(x)(x)"],
            ["", "x", "xxxx"],
            id="recursive-rule-repeat",
        ),
    ],
)
def test_lark_repetition_forms(
    grammar: str, accepted: Sequence[str], rejected: Sequence[str]
) -> None:
    _assert_language(grammar, accepted, rejected)


@pytest.mark.parametrize(
    "common_name, accepted, rejected",
    [
        pytest.param("DIGIT", ["0", "9"], ["", "a", "10"], id="DIGIT"),
        pytest.param("HEXDIGIT", ["0", "a", "F"], ["", "g", "ff"], id="HEXDIGIT"),
        pytest.param("INT", ["0", "123"], ["", "-1", "1.0"], id="INT"),
        pytest.param("SIGNED_INT", ["0", "+12", "-3"], ["", "+", "1.0"], id="SIGNED_INT"),
        pytest.param("DECIMAL", ["1.", "1.5", ".25"], ["", "1", "."], id="DECIMAL"),
        pytest.param("_EXP", ["e1", "E+12", "e-3"], ["", "1", "e"], id="_EXP"),
        pytest.param("FLOAT", ["1.", ".5", "1e3", "1.2e-3"], ["", "1", "e3"], id="FLOAT"),
        pytest.param(
            "SIGNED_FLOAT", ["-1.", "+.5", "-1e3", "1.2e-3"], ["", "1", "+"], id="SIGNED_FLOAT"
        ),
        pytest.param("NUMBER", ["0", "12", ".5", "1e3"], ["", "-1", "x"], id="NUMBER"),
        pytest.param(
            "SIGNED_NUMBER", ["0", "-12", "+.5", "-1e3"], ["", "+", "x"], id="SIGNED_NUMBER"
        ),
        pytest.param(
            "ESCAPED_STRING",
            ['""', '"abc"', '"a\\"b"', '"a\\nb"'],
            ["", "abc", '"unterminated'],
            id="ESCAPED_STRING",
        ),
        pytest.param("LCASE_LETTER", ["a", "z"], ["", "A", "aa"], id="LCASE_LETTER"),
        pytest.param("UCASE_LETTER", ["A", "Z"], ["", "a", "AA"], id="UCASE_LETTER"),
        pytest.param("LETTER", ["a", "Z"], ["", "1", "ab"], id="LETTER"),
        pytest.param("WORD", ["a", "AbCd"], ["", "a1", "a_b"], id="WORD"),
        pytest.param("CNAME", ["a", "_a1", "A_2"], ["", "1a", "a-b"], id="CNAME"),
        pytest.param("WS_INLINE", [" ", "\t", " \t "], ["", "\n", "a"], id="WS_INLINE"),
        pytest.param("WS", [" ", "\n", "\t\f\r\n"], ["", "a"], id="WS"),
        pytest.param("CR", ["\r"], ["", "\n", "\r\n"], id="CR"),
        pytest.param("LF", ["\n"], ["", "\r", "\r\n"], id="LF"),
        pytest.param("NEWLINE", ["\n", "\r\n", "\r\n\n"], ["", "\r", "a"], id="NEWLINE"),
        pytest.param("SH_COMMENT", ["#", "# hello"], ["", "// hello", "# a\n"], id="SH_COMMENT"),
        pytest.param(
            "CPP_COMMENT", ["//", "// hello"], ["", "# hello", "// a\n"], id="CPP_COMMENT"
        ),
        pytest.param(
            "C_COMMENT",
            ["/**/", "/* hello */", "/* ** x **/"],
            ["", "// hello", "/* open"],
            id="C_COMMENT",
        ),
        pytest.param(
            "SQL_COMMENT", ["--", "-- hello"], ["", "- hello", "-- a\n"], id="SQL_COMMENT"
        ),
    ],
)
def test_lark_common_imports(
    common_name: str, accepted: Sequence[str], rejected: Sequence[str]
) -> None:
    _assert_language(f"%import common.{common_name}\nstart: {common_name}", accepted, rejected)


def test_lark_multi_import_alias_and_forward_import() -> None:
    grammar = """
        %ignore WS_INLINE
        start: NAME "=" NUMBER
        %import common (CNAME, WS_INLINE)
        %import common.INT -> NUMBER
        NAME: CNAME
    """
    _assert_language(grammar, ["x=1", "name = 42", "_x\t=\t0"], ["", "1x=2", "x=-1"])


def test_lark_ignore_is_inserted_between_and_after_lexemes() -> None:
    grammar = """
        %import common.WS
        %ignore WS
        start: "a" DIGIT "c".."d"
        DIGIT: "0".."9"
    """
    _assert_language(grammar, ["a1c", "a 1 d", "a\n1\n c  "], [" a1c", "a 1 e", "a x c"])


def test_lark_multiple_ignore_declarations() -> None:
    grammar = """
        %import common (WS, CPP_COMMENT, SH_COMMENT)
        %ignore WS
        %ignore CPP_COMMENT
        %ignore SH_COMMENT
        start: "a" "b"
    """
    _assert_language(
        grammar,
        ["ab", "a // comment\n b", "a# comment\nb", "a b // trailing"],
        [" // initial\na b", "a x b"],
    )


def test_lark_ignore_inline_regex() -> None:
    grammar = r"""
        %ignore /[ _]+/
        start: "a" "b"
    """
    _assert_language(grammar, ["ab", "a_b", "a _ b___"], ["_ab", "a-b"])


def test_lark_allow_initial_skip_options() -> None:
    grammar = """
        %llguidance {"allow_initial_skip": false}
        %llguidance {"allow_initial_skip": true, "no_forcing": false, "allow_invalid_utf8": false}
        %import common.WS
        %ignore WS
        start: "a" "b"
    """
    _assert_language(grammar, ["ab", " a b", "\n\ta\n b  "], ["a c", " xab"])


def test_lark_default_disallows_initial_skip_but_allows_trailing_skip() -> None:
    grammar = """
        %import common.WS_INLINE
        %ignore WS_INLINE
        start: "a" "b"
    """
    _assert_language(grammar, ["ab", "a b", "ab  "], [" ab", "a\nb"])


@pytest.mark.parametrize(
    "schema, accepted, rejected",
    [
        pytest.param(
            '{"type":"string"}',
            ['""', '"hello"', '"λ"'],
            ["hello", "1", '"unterminated'],
            id="string",
        ),
        pytest.param('{"type":"integer"}', ["0", "-12", "123"], ["1.0", '"1"', "+1"], id="integer"),
        pytest.param('{"const":"fixed"}', ['"fixed"'], ['"other"', "fixed"], id="const"),
        pytest.param(
            '{"enum":["red","green",3]}', ['"red"', '"green"', "3"], ['"blue"', "4"], id="enum"
        ),
        pytest.param(
            '{"type":"array","items":{"type":"integer"},"minItems":1,"maxItems":3}',
            ["[1]", "[1,2]", "[ 1, 2, 3 ]"],
            ["[]", "[1,2,3,4]", '["1"]'],
            id="array",
        ),
        pytest.param(
            '{"type":"object","properties":{"x":{"type":"integer"}},"required":["x"],"additionalProperties":false}',
            ['{"x":1}', '{ "x" : -2 }'],
            ["{}", '{"x":"1"}', '{"x":1,"y":2}'],
            id="object",
        ),
        pytest.param(
            '{"anyOf":[{"type":"integer"},{"type":"boolean"}]}',
            ["1", "-2", "true", "false"],
            ['"1"', "null"],
            id="any-of",
        ),
        pytest.param(
            '{"type":"string","pattern":"^a[0-9]+$"}',
            ['"a0"', '"a123"'],
            ['"a"', '"ba1"', '"a1x"'],
            id="string-pattern",
        ),
    ],
)
def test_lark_inline_json_schemas(
    schema: str, accepted: Sequence[str], rejected: Sequence[str]
) -> None:
    _assert_language(f"start: %json {schema}", accepted, rejected)


def test_lark_inline_json_inside_sequence_and_repeat() -> None:
    grammar = r"""
        start: "values=" value (";" value)* "."
        value: %json {"type":"integer"}
    """
    _assert_language(
        grammar, ["values=1.", "values=1;-2;3."], ["values=.", 'values="1".', "values=1;."]
    )


def test_lark_nested_lark_has_an_independent_namespace() -> None:
    grammar = """
        start: item %lark {
          start: item
          item: "b"
        } item
        item: "a"
    """
    _assert_language(grammar, ["aba"], ["aaa", "abb", "bbb"])


def test_lark_nested_lark_supports_recursion_json_and_ignore() -> None:
    grammar = r"""
        start: "[" %lark {
          %llguidance {"allow_initial_skip": true}
          %import common.WS
          %ignore WS
          start: item ":" %json {"type":"integer"}
          item: "x" | "(" item ")"
        } "]"
    """
    _assert_language(grammar, ["[x:1]", "[ ((x)) : -2 ]"], ["[():1]", '[x:"1"]', "[x 1]"])


def test_lark_multiple_nested_grammars() -> None:
    grammar = """
        start: %lark { start: "a" | "b" } %lark { start: "1" | "2" }
    """
    _assert_language(grammar, ["a1", "a2", "b1", "b2"], ["", "a", "1a", "c1"])


def test_lark_numeric_and_named_special_tokens() -> None:
    tokenizer_info = xgr.TokenizerInfo(["a", "<|tool|>", "b", "</s>"], stop_token_ids=[3])
    _assert_token_language(
        "start: <[0,2]> | <|tool|>",
        tokenizer_info,
        accepted=[[0], [1], [2]],
        rejected=[[3], [0, 1], []],
    )


@pytest.mark.parametrize(
    "grammar, accepted, rejected",
    [
        pytest.param("start: <[0-2,1-3,3]>", [[0], [1], [2], [3]], [[4]], id="merged-ranges"),
        pytest.param("start: <[^1,3]>", [[0], [2], [4]], [[1], [3]], id="excluded-set"),
        pytest.param("start: <[*]>", [[0], [1], [2], [3], [4]], [], id="wildcard"),
        pytest.param("start: <[0]> <[2-3]>", [[0, 2], [0, 3]], [[0], [2, 0]], id="token-sequence"),
    ],
)
def test_lark_numeric_special_token_sets(
    grammar: str, accepted: Sequence[Sequence[int]], rejected: Sequence[Sequence[int]]
) -> None:
    tokenizer_info = xgr.TokenizerInfo(["a", "b", "c", "d", "e"])
    _assert_token_language(grammar, tokenizer_info, accepted, rejected)


def test_lark_named_special_token_matches_every_exact_vocab_entry() -> None:
    tokenizer_info = xgr.TokenizerInfo(["<dup>", "x", "<dup>", "dup", "<other>"])
    _assert_token_language(
        "start: <dup>", tokenizer_info, accepted=[[0], [2]], rejected=[[1], [3], [4]]
    )


def test_lark_special_token_and_literal_sequence() -> None:
    tokenizer_info = xgr.TokenizerInfo(["<|tool|>", "x", "y"])
    _assert_token_language(
        'start: <|tool|> "x"', tokenizer_info, accepted=[[0, 1]], rejected=[[0], [0, 2], [1]]
    )


TOOL_CALL_GRAMMAR = r"""
    start: tool* tail
    tail: TEXT

    tool_head[lazy]: TEXT "<tool_call>"
    tool: tool_head %json {
      "type": "object",
      "properties": {"x": {"type": "integer"}},
      "required": ["x"],
      "additionalProperties": false
    } "</tool_call>"

    TEXT: /(\n|.)*/
"""


def test_lark_dynamic_tool_call_optional_repeated_and_committed() -> None:
    _assert_language(
        TOOL_CALL_GRAMMAR,
        [
            "",
            "plain text",
            "text <tool_cal and more",
            '<tool_call>{"x":1}</tool_call>',
            'before<tool_call>{"x":1}</tool_call>after',
            '<tool_call>{"x":1}</tool_call><tool_call>{"x":2}</tool_call>',
            'line 1\nline 2<tool_call>{ "x" : -3 }</tool_call>tail',
        ],
        [
            "<tool_call>",
            '<tool_call>{"x":"bad"}</tool_call>',
            '<tool_call>{"x":1}',
            "before<tool_call>free text</tool_call>after",
            '<tool_call> {"x":1}</tool_call>',
            '<tool_call>{"x":1} </tool_call>',
        ],
    )


def test_lark_dynamic_distinct_string_triggers() -> None:
    grammar = r"""
        start: (foo | bar)* tail
        tail: TEXT

        foo_head[lazy]: TEXT "<foo>"
        foo: foo_head /[a-z]+/ "</foo>"

        bar_head[lazy]: TEXT "<bar>"
        bar: bar_head /[0-9]+/ "</bar>"

        TEXT: /(\n|.)*/
    """
    _assert_language(
        grammar,
        ["free text", "partial <fo remains text", "x<foo>abc</foo>y", "<bar>12</bar><foo>x</foo>"],
        ["<foo>", "<foo>12</foo>", "<bar>x</bar>", "<bar>12"],
    )


def test_lark_dynamic_shared_trigger_dispatch() -> None:
    grammar = r"""
        start: (foo | bar)* tail
        tail: TEXT

        foo_head[lazy]: TEXT "<function"
        foo: foo_head "=foo>" /[a-z]+/ "</function>"

        bar_head[lazy]: TEXT "<function"
        bar: bar_head "=bar>" /[A-Z]+/ "</function>"

        TEXT: /(\n|.)*/
    """
    _assert_language(
        grammar,
        [
            "free text",
            "a<function=foo>abc</function>b",
            "<function=bar>ABC</function><function=foo>xyz</function>",
        ],
        [
            "<function",
            "<function=baz>abc</function>",
            "<function=foo>ABC</function>",
            "<function=bar>abc</function>",
            "<function=foo>abc",
        ],
    )


def test_lark_dynamic_any_text_can_be_referenced_through_terminals() -> None:
    grammar = r"""
        start: tool* tail
        tail: FREE
        head[lazy]: FREE "<call>"
        tool: head /[0-9]+/ "</call>"
        FREE: TEXT
        TEXT: /(\n|.)*/
    """
    _assert_language(
        grammar,
        ["free", "x<call>12</call>y", "partial <cal"],
        ["<call>", "<call>x</call>", "<call>12"],
    )


def test_lark_standalone_lazy_rule() -> None:
    grammar = r"""
        start: head
        head[lazy]: TEXT "<end>"
        TEXT: /(\n|.)*/
    """
    _assert_language(grammar, ["", "plain", "<end>", "plain<end>"], ["<end>x", "a<end>b"])


def test_lark_dynamic_special_token_trigger() -> None:
    tokenizer_info = xgr.TokenizerInfo(
        ["plain", "<|tool|>", "{", '"x"', ":", "1", "}", "</tool>", "bad", "</s>"],
        stop_token_ids=[9],
    )
    grammar = r"""
        start: tool* tail
        tail: TEXT
        tool: TEXT <|tool|> %json {
          "type": "object",
          "properties": {"x": {"const": 1}},
          "required": ["x"],
          "additionalProperties": false
        } "</tool>"
        TEXT: /(\n|.)*/
    """
    _assert_token_language(
        grammar,
        tokenizer_info,
        accepted=[[0], [1, 2, 3, 4, 5, 6, 7], [0, 1, 2, 3, 4, 5, 6, 7, 0]],
        rejected=[[1], [1, 8], [1, 2, 3, 4, 8]],
    )


def test_lark_serialization_round_trip_for_core_and_dynamic_grammars() -> None:
    core = xgr.Grammar.from_lark('start: "a" ("b" | "c")?')
    restored_core = xgr.Grammar.deserialize_json(core.serialize_json())
    _assert_grammar_language(restored_core, ["a", "ab", "ac"], ["", "abc"])

    dynamic = xgr.Grammar.from_lark(TOOL_CALL_GRAMMAR)
    restored_dynamic = xgr.Grammar.deserialize_json(dynamic.serialize_json())
    _assert_grammar_language(
        restored_dynamic,
        ["text", '<tool_call>{"x":1}</tool_call>tail'],
        ["<tool_call>", '<tool_call>{"x":"bad"}</tool_call>'],
    )


def test_lark_serialization_round_trip_for_token_dispatch() -> None:
    tokenizer_info = xgr.TokenizerInfo(["plain", "<|call|>", "x", "</call>"])
    grammar = xgr.Grammar.from_lark(
        r"""
        start: call* tail
        tail: TEXT
        call: TEXT <|call|> "x" "</call>"
        TEXT: /(\n|.)*/
        """,
        tokenizer_info=tokenizer_info,
    )
    restored = xgr.Grammar.deserialize_json(grammar.serialize_json())
    compiled = xgr.GrammarCompiler(tokenizer_info, cache_enabled=False).compile_grammar(restored)
    assert _matches_token_sequence(compiled, [0])
    assert _matches_token_sequence(compiled, [1, 2, 3])
    assert not _matches_token_sequence(compiled, [1, 3])


def test_lark_grammar_union_and_concat_integration() -> None:
    left = xgr.Grammar.from_lark('start: "a" | "b"')
    right = xgr.Grammar.from_lark("start: /[0-9]+/")

    _assert_grammar_language(xgr.Grammar.union(left, right), ["a", "b", "0", "123"], ["a1", "c"])
    _assert_grammar_language(xgr.Grammar.concat(left, right), ["a0", "b123"], ["a", "1", "c1"])


def test_lark_large_choice_grammar() -> None:
    options = [f"option-{index:03d}" for index in range(256)]
    choices = " | ".join(json.dumps(option) for option in options)
    grammar = f"start: option\noption: {choices}"
    _assert_language(
        grammar, [options[0], options[127], options[-1]], ["", "option-256", "option-12"]
    )


@pytest.mark.parametrize(
    "grammar, message",
    [
        pytest.param('item: "a"', "no start rule", id="missing-start"),
        pytest.param("start: missing", "unknown name 'missing'", id="unknown-rule"),
        pytest.param(
            'start: foo\nfoo: "a"\nfoo: "b"',
            "duplicate rule or terminal 'foo'",
            id="duplicate-rule",
        ),
        pytest.param(
            'start: FOO\nFOO: "a"\nFOO: "b"',
            "duplicate rule or terminal 'FOO'",
            id="duplicate-terminal",
        ),
        pytest.param(
            "start: FOO\nFOO: BAR\nBAR: FOO", "circular reference in terminal", id="terminal-cycle"
        ),
        pytest.param(
            'start: TOKEN\nTOKEN: rule\nrule: "a"',
            "terminal 'TOKEN' cannot reference rule 'rule'",
            id="terminal-references-rule",
        ),
        pytest.param(
            'start[capture]: "a"', "attribute 'capture' is not supported", id="capture-attribute"
        ),
        pytest.param(
            'start[stop="x"]: "a"', "attribute 'stop' is not supported", id="stop-attribute"
        ),
        pytest.param(
            "TOKEN[lazy]: /a/\nstart: TOKEN",
            "attributes are only supported on rules",
            id="terminal-attribute",
        ),
        pytest.param('start.1: "a"', "priorities are not supported", id="priority"),
        pytest.param('start{x}: "a"', "Lark templates are not supported", id="template-definition"),
        pytest.param(
            'start: foo{x}\nfoo: "a"', "Lark templates are not supported", id="template-reference"
        ),
        pytest.param(
            'start::0: "a"', "parametric grammar is not supported", id="parametric-definition"
        ),
        pytest.param(
            "start: foo::0", "parametric grammar is not supported", id="parametric-reference"
        ),
        pytest.param(
            'start: "a" %if enabled',
            "parametric %if conditions are not supported",
            id="parametric-if",
        ),
        pytest.param(
            'start: A & B\nA: "a"\nB: "b"', "intersection '&' is not supported", id="intersection"
        ),
        pytest.param("start: ~/[a]/", "complement '~' is not supported", id="complement"),
        pytest.param(
            'start: "a"i', "case-insensitive string flags are not supported", id="string-flag"
        ),
        pytest.param(
            "start: /abc/i", "regular-expression flags are not supported", id="regex-flag"
        ),
        pytest.param("start: /[abc/", "failed to compile regular expression", id="invalid-regex"),
        pytest.param('start: "\\q"', "invalid string literal", id="invalid-string-escape"),
        pytest.param(
            'start: "unterminated', "unterminated string literal", id="unterminated-string"
        ),
        pytest.param(
            "start: /unterminated", "unterminated regular expression", id="unterminated-regex"
        ),
        pytest.param(
            "start: <unterminated", "unterminated special token", id="unterminated-special-token"
        ),
        pytest.param("start: <bad token>", "invalid special token", id="special-token-whitespace"),
        pytest.param("start: @", "empty grammar reference", id="empty-grammar-reference"),
        pytest.param("start: $", "unexpected character '$'", id="unexpected-character"),
        pytest.param("start: -", "unexpected '-' character", id="unexpected-minus"),
        pytest.param("start: !", "unexpected '!' character", id="unexpected-bang"),
        pytest.param('start "a"', "expected ':' after rule name", id="missing-colon"),
        pytest.param('start: ("a"', "expected ')' after group", id="unclosed-group"),
        pytest.param('start: ["a"', "expected ']' after optional group", id="unclosed-optional"),
        pytest.param('start: "a" ->', "expected alias name after '->'", id="missing-alias-name"),
        pytest.param(
            "%declare TOKEN\nstart: TOKEN",
            "directive %declare is not supported",
            id="declare-directive",
        ),
        pytest.param(
            '%override start\nstart: "a"',
            "directive %override is not supported",
            id="override-directive",
        ),
        pytest.param(
            "%import common.UNKNOWN\nstart: UNKNOWN",
            "unknown common import",
            id="unknown-common-import",
        ),
        pytest.param(
            '%import common.INT\nINT: "x"\nstart: INT',
            "duplicate rule or terminal 'INT'",
            id="import-name-conflict",
        ),
        pytest.param('start: "a"{3,2}', "repetition end must be greater", id="repetition-reversed"),
        pytest.param('start: "a"{-1,}', "invalid repetition count", id="negative-repetition"),
        pytest.param(
            'start: "a"{999999999999999999999}',
            "invalid repetition count",
            id="repetition-overflow",
        ),
        pytest.param(
            'start: "ab".."c"', "range endpoints must be one character", id="range-start-too-long"
        ),
        pytest.param(
            'start: "a".."bc"', "range endpoints must be one character", id="range-end-too-long"
        ),
        pytest.param('start: "z".."a"', "range start must not exceed end", id="range-reversed"),
        pytest.param(
            "start: %json {",
            "failed to parse JSON value after %json",
            id="malformed-json-directive",
        ),
        pytest.param(
            "start: %json []", "failed to compile inline JSON schema", id="invalid-json-schema"
        ),
        pytest.param(
            'start: %lark { item: "a" }',
            "failed to compile nested Lark grammar",
            id="nested-no-start",
        ),
        pytest.param(
            'start: %regex {"substring_chars":"abc"}',
            "structured %regex is not supported",
            id="structured-regex",
        ),
        pytest.param(
            "start: @other",
            "multiple grammar references are not supported",
            id="named-grammar-reference",
        ),
        pytest.param(
            "start: TOKEN\nTOKEN: <[1]>",
            "special tokens cannot be used in terminals",
            id="special-in-terminal",
        ),
        pytest.param(
            "start: <[1-2-3]>", "invalid numeric special-token range", id="multiple-range-dashes"
        ),
        pytest.param(
            "start: <[3-1]>", "invalid numeric special-token range", id="numeric-range-reversed"
        ),
        pytest.param("start: <[,]>", "empty numeric special-token range", id="empty-token-range"),
        pytest.param(
            "start: <[*]>",
            "wildcard special token requires tokenizer_info",
            id="wildcard-needs-tokenizer",
        ),
        pytest.param(
            "start: <[^*]>",
            "negated wildcard special token is not supported",
            id="negated-wildcard",
        ),
        pytest.param("start: <[*,1]>", "wildcard cannot be mixed", id="mixed-wildcard-range"),
        pytest.param(
            "start: <|tool|>",
            "named special token <|tool|> requires tokenizer_info",
            id="named-needs-tokenizer",
        ),
        pytest.param(
            "start: <[0-1000001]>", "special-token range is too large", id="token-range-too-large"
        ),
        pytest.param(
            '%llguidance {"allow_initial_skip": 1}\nstart: "a"',
            "allow_initial_skip must be a boolean",
            id="initial-skip-type",
        ),
        pytest.param(
            '%llguidance {"no_forcing": true}\nstart: "a"',
            "%llguidance option 'no_forcing' is not supported",
            id="no-forcing-option",
        ),
        pytest.param(
            '%llguidance {"allow_invalid_utf8": true}\nstart: "a"',
            "%llguidance option 'allow_invalid_utf8' is not supported",
            id="invalid-utf8-option",
        ),
        pytest.param(
            '%llguidance {"unknown": false}\nstart: "a"',
            "unknown %llguidance option 'unknown'",
            id="unknown-llguidance-option",
        ),
        pytest.param(
            '%llguidance []\nstart: "a"',
            "%llguidance value must be an object",
            id="llguidance-not-object",
        ),
        pytest.param(
            "start: thing\nthing[lazy]: /[a-z]+/",
            "general lazy rules are not supported",
            id="unsupported-general-lazy",
        ),
        pytest.param(
            '%ignore MISSING\nstart: "a"', "unknown name 'MISSING'", id="unknown-ignore-name"
        ),
    ],
)
def test_lark_errors_are_explicit_and_located(grammar: str, message: str) -> None:
    _assert_lark_error(grammar, message)


def test_lark_named_special_token_error_with_tokenizer() -> None:
    tokenizer_info = xgr.TokenizerInfo(["<known>", "text"])
    _assert_lark_error("start: <unknown>", "unknown special token <unknown>", tokenizer_info)


def test_lark_dynamic_trigger_levels_cannot_be_mixed() -> None:
    tokenizer_info = xgr.TokenizerInfo(["<|bar|>", "x"])
    grammar = r"""
        start: (foo | bar)* tail
        tail: TEXT
        foo_head[lazy]: TEXT "<foo>"
        foo: foo_head "x"
        bar: TEXT <|bar|> "x"
        TEXT: /(\n|.)*/
    """
    _assert_lark_error(grammar, "cannot mix string and token triggers", tokenizer_info)


def test_lark_lazy_and_dynamic_special_token_triggers_cannot_be_negated() -> None:
    tokenizer_info = xgr.TokenizerInfo(["a", "b", "c"])
    _assert_lark_error(
        "start: head\nhead[lazy]: TEXT <[^1]>\nTEXT: /(\\n|.)*/",
        "lazy special-token trigger cannot be negated",
        tokenizer_info,
    )
    _assert_lark_error(
        'start: tool* tail\ntail: TEXT\ntool: TEXT <[^1]> "x"\nTEXT: /(\\n|.)*/',
        "dynamic special-token trigger cannot be negated",
        tokenizer_info,
    )


def test_lark_error_reports_crlf_line_column_and_source_context() -> None:
    error = _assert_lark_error(
        '# comment\r\nstart: "a"\r\nitem missing', "expected ':' after rule name"
    )
    assert "line 3, column 6" in error
    assert "item missing" in error
    assert "     ^" in error


@pytest.mark.parametrize(
    "grammar, accepted, rejected",
    [
        pytest.param('start: ("a"?)*', ["", "a", "aaaa"], ["b", "ab"], id="nullable-star"),
        pytest.param(
            'start: ("a"?){2,3}', ["", "a", "aa", "aaa"], ["aaaa", "b"], id="nullable-brace-repeat"
        ),
        pytest.param(
            'start: %json {"type":"boolean"} {2}',
            ["truefalse", "falsefalse"],
            ["true", "truefalsetrue"],
            id="repetition-of-inline-json",
        ),
        pytest.param(
            'start: item {2}\nitem: "a"', ["aa"], ["a", "aaa"], id="brace-repetition-after-space"
        ),
        pytest.param('?start: "a"', ["a"], ["", "aa"], id="question-prefixed-start"),
        pytest.param('!start: "a"', ["a"], ["", "aa"], id="bang-prefixed-start"),
        pytest.param(
            'start: ("a" | "b") -> outer', ["a", "b"], ["outer", "ab"], id="alias-after-group"
        ),
        pytest.param(
            'start: ( "a"\n  | "b"\n  )', ["a", "b"], ["", "ab"], id="multiline-group-body"
        ),
        pytest.param(
            'start: [ "a"\n ] "b"', ["ab", "b"], ["", "a"], id="multiline-optional-group-body"
        ),
        pytest.param('start: "a"\r\n  | "b"\r\n', ["a", "b"], ["", "ab"], id="crlf-line-endings"),
        pytest.param(r'start: "\/" "A"', ["/A"], ["\\/A", "/"], id="escaped-slash"),
        pytest.param(
            'start: "\\ud83d\\ude00"',
            ["\U0001f600"],
            ["", "\\ud83d\\ude00"],
            id="surrogate-pair-escape",
        ),
        pytest.param(
            'start: "😀".."😂"',
            ["😀", "😁", "😂"],
            ["", "😃", "a", "😀😁"],
            id="supplementary-plane-range",
        ),
        pytest.param(
            "start: NUM\n%import common.DIGIT\nNUM: DIGIT DIGIT",
            ["12", "00"],
            ["", "1", "123"],
            id="terminal-composed-of-imported-terminal",
        ),
    ],
)
def test_lark_additional_core_languages(
    grammar: str, accepted: Sequence[str], rejected: Sequence[str]
) -> None:
    _assert_language(grammar, accepted, rejected)


def test_lark_deeply_nested_groups() -> None:
    depth = 400
    grammar = "start: " + "(" * depth + '"a"' + ")" * depth
    _assert_language(grammar, ["a"], ["", "aa"])


def test_lark_mutually_empty_recursion_matches_nothing() -> None:
    grammar = "start: a\na: b\nb: a"
    _assert_language(grammar, [], ["", "a", "x"])


def test_lark_self_recursive_start_matches_nothing() -> None:
    _assert_language("start: start", [], ["", "a"])


@pytest.mark.parametrize(
    "grammar, accepted, rejected",
    [
        pytest.param(
            """
            %import common.WS
            %ignore WS
            start: WS "a"
            """,
            [" a", "  a", "\ta "],
            ["a", "b"],
            id="ignored-terminal-also-referenced-directly",
        ),
        pytest.param(
            """
            %ignore PAD
            PAD: ("_" | "-"){1,2}
            start: "a" "b"
            """,
            ["ab", "a_b", "a-_b", "a___b"],
            ["_ab", "a b"],
            id="composed-ignore-terminal-loops",
        ),
        pytest.param(
            """
            %import common.WS
            %ignore WS
            start: "a".."z" "0".."9"
            """,
            ["a1", "a 1", "a1 "],
            [" a1", "a1x"],
            id="character-ranges-are-lexemes-for-skip",
        ),
        pytest.param(
            """
            %import common.WS
            %ignore WS
            start: ITEM{2}
            ITEM: "x"
            """,
            ["xx", "x x", "x x "],
            ["x", " xx"],
            id="repeated-terminal-allows-skip-between",
        ),
        pytest.param(
            """
            %import common.WS -> W
            %ignore W
            start: "a" "b"
            """,
            ["ab", "a b"],
            [" ab", "a-b"],
            id="ignore-through-import-alias",
        ),
    ],
)
def test_lark_ignore_edge_cases(
    grammar: str, accepted: Sequence[str], rejected: Sequence[str]
) -> None:
    _assert_language(grammar, accepted, rejected)


def test_lark_dynamic_prefix_overlapping_triggers_are_both_reachable() -> None:
    grammar = r"""
        start: (foo | bar)* tail
        tail: TEXT

        foo_head[lazy]: TEXT "<a"
        foo: foo_head ">" /[0-9]+/ "</>"

        bar_head[lazy]: TEXT "<ab"
        bar: bar_head ">" /[x-z]+/ "</>"

        TEXT: /(\n|.)*/
    """
    _assert_language(
        grammar,
        ["free", "<a>12</>", "<ab>x</>", "pre<a>1</>mid<ab>z</>post"],
        ["<ab>12</>", "<a>x</>", "<a>"],
    )


def test_lark_dynamic_unicode_trigger_strings() -> None:
    grammar = 'start: tool* tail\ntail: TEXT\nhead[lazy]: TEXT "工具>"\n'
    grammar += 'tool: head /[0-9]+/ "</工具>"\nTEXT: /(\\n|.)*/'
    _assert_language(
        grammar, ["文本工具>42</工具>后面", "工具>1</工具>"], ["工具>x</工具>", "工具>42"]
    )


def test_lark_dynamic_tool_body_may_reference_the_tail_rule() -> None:
    # The rules consumed by the dynamic start lowering keep their real bodies, so
    # a tool body may reference them and still match the written language.
    grammar = r"""
        start: tool* tail
        tail: TEXT
        head[lazy]: TEXT "<t>"
        tool: head "A" tail "B"
        TEXT: /(\n|.)*/
    """
    _assert_language(
        grammar, ["x<t>ABy", "x<t>A hello By", "<t>A\nmulti\nline\nB"], ["x<t>A", "x<t>B"]
    )


def test_lark_dynamic_tool_body_may_reference_tools_recursively() -> None:
    grammar = r"""
        start: tool* tail
        tail: TEXT
        head[lazy]: TEXT "<t>"
        tool: head "A" inner "B"
        inner: tool | "i"
        TEXT: /(\n|.)*/
    """
    _assert_language(grammar, ["x<t>AiBy", "x<t>A<t>AiBBy"], ["x<t>ABy", "x<t>A<t>ABBy"])


def test_lark_dynamic_allows_ignored_content_after_string_trigger() -> None:
    grammar = r"""
        %import common.WS
        %ignore WS
        start: tool* tail
        tail: TEXT
        head[lazy]: TEXT "<t>"
        tool: head "a" "b" "</t>"
        TEXT: /(\n|.)*/
    """
    _assert_language(
        grammar,
        ["x<t>ab</t>y", "x<t>a b </t>y", "x<t> ab</t>y", "x<t>\na\nb\n</t>y"],
        ["x<t>axb</t>y", "x<t>ab"],
    )


def test_lark_dynamic_allows_ignored_content_after_token_trigger() -> None:
    tokenizer_info = xgr.TokenizerInfo(["<|call|>", "x", " ", "</c>", "plain"])
    grammar = r"""
        %ignore / +/
        start: call* tail
        tail: TEXT
        call: TEXT <|call|> "x" "</c>"
        TEXT: /(\n|.)*/
    """
    _assert_token_language(
        grammar,
        tokenizer_info,
        accepted=[[4], [0, 1, 3], [0, 2, 1, 3], [0, 1, 2, 3]],
        rejected=[[0, 3], [0, 1], [0, 4, 1, 3]],
    )


def test_lark_dynamic_tail_can_be_a_terminal_reference() -> None:
    grammar = r"""
        start: tool* TEXT
        head[lazy]: TEXT "<c>"
        tool: head /[0-9]+/ "</c>"
        TEXT: /(\n|.)*/
    """
    _assert_language(grammar, ["free", "x<c>12</c>y", "<c>1</c><c>2</c>"], ["<c>", "<c>x</c>"])


def test_lark_dynamic_tool_with_empty_remainder() -> None:
    grammar = r"""
        start: tool* tail
        tail: TEXT
        head[lazy]: TEXT "<t>"
        tool: head
        TEXT: /(\n|.)*/
    """
    _assert_language(grammar, ["", "plain", "x<t>y<t>z", "<t>"], [])


def test_lark_non_dynamic_start_with_extra_element_falls_back_to_lazy_rules() -> None:
    grammar = r"""
        start: tool* tail "!"
        tail: TEXT
        head[lazy]: TEXT "<t>"
        tool: head /[0-9]+/ "</t>"
        TEXT: /(\n|.)*/
    """
    _assert_language(grammar, ["plain!", "x<t>1</t>y!", "x<t>1</t>!"], ["plain", "x<t>1"])


def test_lark_inline_dotall_group_regex_works_as_lazy_text() -> None:
    grammar = 'start: head\nhead[lazy]: /(?s:.*)/ "<end>"'
    _assert_language(grammar, ["ab", "ab<end>", "<end>"], ["<end>x", "a<end>b"])


def test_lark_lazy_start_rule_itself() -> None:
    grammar = 'start[lazy]: TEXT "<end>"\nTEXT: /(\\n|.)*/'
    _assert_language(grammar, ["", "plain", "<end>", "plain<end>"], ["<end>x", "a<end>b"])


def test_lark_lazy_rule_commits_at_the_first_trigger() -> None:
    grammar = 'start: head "x"\nhead[lazy]: TEXT "<end>"\nTEXT: /(\\n|.)*/'
    _assert_language(grammar, ["abc<end>x", "<end>x"], ["abc<end>", "a<end>b<end>x"])


def test_lark_dead_rules_are_pruned_from_dynamic_grammars() -> None:
    grammar = xgr.Grammar.from_lark(
        r"""
        start: tool* tail
        tail: TEXT
        tool_head[lazy]: TEXT "<tool_call>"
        tool: tool_head /[0-9]+/ "</tool_call>"
        TEXT: /(\n|.)*/
        """
    )
    printed = str(grammar)
    assert "TagDispatch" in printed
    assert "tool_head" not in printed
    assert "tail ::=" not in printed
    assert "TEXT ::=" not in printed


def test_lark_out_of_vocab_token_ids_never_match_and_do_not_corrupt_masks() -> None:
    # The grammar is converted without tokenizer metadata, then compiled against a
    # small vocab. Out-of-vocab ids must be ignored rather than corrupt the mask.
    tokenizer_info = xgr.TokenizerInfo(["a", "b", "c", "d"])
    grammar = xgr.Grammar.from_lark("start: <[2,100]>")
    compiled = xgr.GrammarCompiler(tokenizer_info, cache_enabled=False).compile_grammar(grammar)

    matcher = xgr.GrammarMatcher(compiled, terminate_without_stop_token=True)
    bitmask = xgr.allocate_token_bitmask(1, tokenizer_info.vocab_size)
    matcher.fill_next_token_bitmask(bitmask)
    rejected = set(_get_masked_tokens_from_bitmask(bitmask, tokenizer_info.vocab_size))
    assert set(range(tokenizer_info.vocab_size)) - rejected == {2}

    assert _matches_token_sequence(compiled, [2])
    assert not _matches_token_sequence(compiled, [0])


def test_lark_int32_max_token_id_compiles_against_small_vocab() -> None:
    tokenizer_info = xgr.TokenizerInfo(["a", "b"])
    grammar = xgr.Grammar.from_lark("start: <[2147483647]>")
    compiled = xgr.GrammarCompiler(tokenizer_info, cache_enabled=False).compile_grammar(grammar)
    assert not _matches_token_sequence(compiled, [0])
    assert not _matches_token_sequence(compiled, [1])


def test_lark_numeric_token_ids_are_validated_against_tokenizer_when_provided() -> None:
    tokenizer_info = xgr.TokenizerInfo(["a", "b"])
    _assert_lark_error(
        "start: <[100]>", "out of range for vocab size 2", tokenizer_info=tokenizer_info
    )
    _assert_lark_error(
        "start: <[0-100]>", "out of range for vocab size 2", tokenizer_info=tokenizer_info
    )
    _assert_lark_error(
        "start: <[^0,100]>", "out of range for vocab size 2", tokenizer_info=tokenizer_info
    )


def test_lark_numeric_token_single_range_and_leading_zeros() -> None:
    tokenizer_info = xgr.TokenizerInfo(["a", "b", "c"])
    _assert_token_language(
        "start: <[0-0]> <[01]>", tokenizer_info, accepted=[[0, 1]], rejected=[[0, 0], [1, 1]]
    )


def test_lark_special_token_sequence_allows_ignored_content() -> None:
    tokenizer_info = xgr.TokenizerInfo(["<|t|>", "x", " "])
    _assert_token_language(
        '%ignore / +/\nstart: <|t|> "x"',
        tokenizer_info,
        accepted=[[0, 1], [0, 2, 1]],
        rejected=[[1], [0, 2]],
    )


def test_lark_excluded_token_set_serialization_round_trip() -> None:
    tokenizer_info = xgr.TokenizerInfo(["a", "b", "c", "d"])
    grammar = xgr.Grammar.from_lark("start: <[^1]>")
    restored = xgr.Grammar.deserialize_json(grammar.serialize_json())
    compiled = xgr.GrammarCompiler(tokenizer_info, cache_enabled=False).compile_grammar(restored)

    matcher = xgr.GrammarMatcher(compiled, terminate_without_stop_token=True)
    bitmask = xgr.allocate_token_bitmask(1, tokenizer_info.vocab_size)
    matcher.fill_next_token_bitmask(bitmask)
    rejected = set(_get_masked_tokens_from_bitmask(bitmask, tokenizer_info.vocab_size))
    assert set(range(tokenizer_info.vocab_size)) - rejected == {0, 2, 3}

    assert _matches_token_sequence(compiled, [0])
    assert not _matches_token_sequence(compiled, [1])


def test_lark_token_level_grammars_do_not_accept_strings() -> None:
    tokenizer_info = xgr.TokenizerInfo(["a", "b"])
    compiled = _compile_lark("start: <[0]>", tokenizer_info)
    matcher = xgr.GrammarMatcher(compiled, terminate_without_stop_token=True)
    assert not matcher.accept_string("a")


@pytest.mark.parametrize(
    "grammar, message",
    [
        pytest.param("", "no start rule found", id="empty-grammar"),
        pytest.param("# only a comment\n", "no start rule found", id="comment-only-grammar"),
        pytest.param('start: "a" ~ b\nb: "x"', "expected integer", id="tilde-followed-by-name"),
        pytest.param('start: "a"{1.5}', "invalid repetition count", id="fractional-repetition"),
        pytest.param('start: "a"~2..', "expected integer", id="tilde-range-missing-end"),
        pytest.param('start: "a\nb"', "unterminated string literal", id="newline-inside-string"),
        pytest.param(r'start: "\x41"', "invalid string literal", id="hex-escape-not-json"),
        pytest.param(
            '%ignore <[1]>\nstart: "a"',
            "special tokens cannot be used in terminals",
            id="ignore-of-special-token",
        ),
        pytest.param(
            'start: "a" item: "b"', "expected grammar expression", id="two-definitions-one-line"
        ),
        pytest.param(
            'start: "a"\n  "b"', "expected rule or terminal name", id="continuation-without-pipe"
        ),
        pytest.param(
            'start: ("a" |\n "b")', "expected ')' after group", id="pipe-at-end-of-line-in-group"
        ),
        pytest.param(
            'start: "a"i.."c"',
            "case-insensitive string flags are not supported",
            id="range-start-flag",
        ),
        pytest.param(
            'start: "a".."c"i',
            "case-insensitive string flags are not supported",
            id="range-end-flag",
        ),
        pytest.param(
            'start: head\nhead[lazy]: TEXT ""\nTEXT: /(\\n|.)*/',
            "general lazy rules are not supported",
            id="empty-string-lazy-trigger",
        ),
        pytest.param(
            'start: head\nhead[lazy]: TEXT %json {"type":"integer"}\nTEXT: /(\\n|.)*/',
            "general lazy rules are not supported",
            id="json-lazy-trigger",
        ),
        pytest.param(
            'start: head\nhead[lazy]: TEXT "<end>"\nTEXT: /(. |\\n)*/',
            "general lazy rules are not supported",
            id="whitespace-in-regex-is-not-any-text",
        ),
        pytest.param(
            '%import common.INT -> item\nitem: "a"\nstart: item',
            "duplicate rule or terminal 'item'",
            id="import-alias-collides-with-rule",
        ),
        pytest.param(
            "%import common.sub.INT\nstart: INT",
            "unknown common import",
            id="nested-unknown-import-path",
        ),
        pytest.param(
            "start: TEXT\nTEXT: /(?s:.*)/",
            "failed to compile regular expression",
            id="dotall-group-terminal-is-rejected-by-regex-converter",
        ),
    ],
)
def test_lark_additional_errors(grammar: str, message: str) -> None:
    _assert_lark_error(grammar, message)


@pytest.mark.parametrize(
    "grammar, alphabet, max_len",
    [
        pytest.param('start: "a" "b" | "c" ("d" | "e")', "abcde", 4, id="choice-groups"),
        pytest.param('start: "x" | "(" start ")"', "x()", 5, id="recursive-start"),
        pytest.param('start: ("ab" | "c")~2..4', "abc", 6, id="tilde-repetition"),
        pytest.param('start: ("a".."c" | "z")+', "abcz", 4, id="ranges-plus"),
        pytest.param('start: ("a"?)* "b"', "ab", 4, id="nullable-star-then-literal"),
        pytest.param("start: A B\nA: /a+/\nB: /ab/", "ab", 5, id="overlapping-terminals"),
        pytest.param("start: /a+/ /ab/ /b+/", "ab", 6, id="ambiguous-regex-split"),
        pytest.param('start: item ("," item)*\nitem: /[ab]+/', "ab,", 5, id="separated-list"),
        pytest.param('start: "a" [ "b" ] "c"', "abc", 4, id="optional-group"),
        pytest.param('start: (A | B)+\nA: "aa"\nB: "ab"', "ab", 6, id="terminal-choice-plus"),
        pytest.param("start: /a|b/ /(ab)+/", "ab", 5, id="alternation-then-group-plus"),
        pytest.param('start: "a"~2..3 "b"?', "ab", 5, id="tilde-range-then-optional"),
        pytest.param(
            'start: e\ne: "a" | "(" f ")"\nf: e ("+" e)*', "a(+)", 5, id="expression-grammar"
        ),
        pytest.param("start: /[ab]*a[ab]/", "ab", 5, id="nondeterministic-regex"),
    ],
)
def test_lark_language_matches_reference_lark_implementation(
    grammar: str, alphabet: str, max_len: int
) -> None:
    """Differential test against the reference `lark` package.

    Enumerates every string over a small alphabet up to a length bound and requires
    xgrammar and lark (earley with the scannerless dynamic_complete lexer) to agree on
    membership. Grammars avoid %ignore, whose placement semantics intentionally follow
    LLGuidance rather than lark.
    """
    lark = pytest.importorskip("lark")

    reference = lark.Lark(grammar, start="start", parser="earley", lexer="dynamic_complete")
    compiled = _compile_lark(grammar)

    for length in range(max_len + 1):
        for candidate in itertools.product(alphabet, repeat=length):
            value = "".join(candidate)
            try:
                reference.parse(value)
                expected = True
            except lark.exceptions.LarkError:
                expected = False
            assert _matches_string(compiled, value) == expected, value
