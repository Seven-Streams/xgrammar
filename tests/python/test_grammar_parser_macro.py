"""Tests the macro features of the grammar parser."""

import string
import sys
import time
from typing import Optional

import pytest
from transformers import AutoConfig, AutoTokenizer

import xgrammar as xgr
from xgrammar.testing import (
    GrammarFunctor,
    _ebnf_to_grammar_no_normalization,
    _is_grammar_accept_string,
)


def test_tag_dispatch():
    """Test TagDispatch functionality."""
    before = """root ::= TagDispatch(
    ("tag1", rule1),
    ("tag2", rule2),
    excludes = ("abc", "def"),
    loop_after_dispatch = false
)
rule1 ::= "a"
rule2 ::= "b"
"""

    expected = """root ::= ((TagDispatch(
  ("tag1", rule1),
  ("tag2", rule2),
  loop_after_dispatch=false,
  excludes=("abc", "def")
)))
rule1 ::= (("a"))
rule2 ::= (("b"))
"""
    grammar = _ebnf_to_grammar_no_normalization(before)
    after = str(grammar)
    assert after == expected


def test_tag_dispatch_default_parameters():
    """Test TagDispatch functionality."""
    before = """root ::= TagDispatch(("tag1", rule1), ("tag2", rule2))
rule1 ::= "a"
rule2 ::= "b"
"""
    expected = """root ::= ((TagDispatch(
  ("tag1", rule1),
  ("tag2", rule2),
  loop_after_dispatch=true,
  excludes=()
)))
rule1 ::= (("a"))
rule2 ::= (("b"))
"""
    grammar = _ebnf_to_grammar_no_normalization(before)
    after = str(grammar)
    assert after == expected


def test_lookahead_assertion_analyzer_tag_dispatch():
    # tag dispatch disables lookahead assertion detection
    before = r"""root ::= TagDispatch(("tag1", rule1), ("tag2", rule2), ("tag3", rule3), ("tag4", rule4), ("tag5", rule5))
rule1 ::= "b"
rule2 ::= "c"
rule3 ::= "" | "d" rule3
rule4 ::= "" | "e" rule4 "f"
rule5 ::= "" | "g" rule5 "h"
"""
    expected = r"""root ::= TagDispatch(
  ("tag1", rule1),
  ("tag2", rule2),
  ("tag3", rule3),
  ("tag4", rule4),
  ("tag5", rule5),
  loop_after_dispatch=true,
  excludes=()
)
rule1 ::= (("b"))
rule2 ::= (("c"))
rule3 ::= ("" | ("d" rule3))
rule4 ::= ("" | ("e" rule4 "f"))
rule5 ::= ("" | ("g" rule5 "h"))
"""

    grammar = _ebnf_to_grammar_no_normalization(before)
    grammar = GrammarFunctor.structure_normalizer(grammar)
    grammar = GrammarFunctor.byte_string_fuser(grammar)
    grammar = GrammarFunctor.lookahead_assertion_analyzer(grammar)
    after = str(grammar)
    assert after == expected


def test_tag_dispatch_end_to_end():
    before = """root ::= TagDispatch(("tag1", rule1), ("tag2", rule2), ("tag3", rule3))
rule1 ::= "a"
rule2 ::= "b"
rule3 ::= "c"
"""
    expected = """root ::= TagDispatch(
  ("tag1", rule1),
  ("tag2", rule2),
  ("tag3", rule3),
  loop_after_dispatch=true,
  excludes=()
)
rule1 ::= (("a"))
rule2 ::= (("b"))
rule3 ::= (("c"))
"""
    grammar = xgr.Grammar.from_ebnf(before)
    after = str(grammar)
    assert after == expected


def test_tag_dispatch_end_to_end_complex():
    before = """root ::= TagDispatch(("tag1", rule1), ("tag2", rule2), ("tag3", rule3))
rule1 ::= ("a" TagDispatch(("tag1", rule2), ("tag2", rule3)) | "zzz")
rule2 ::= TagDispatch(("tag1", rule2), ("tag2", rule3)) | TagDispatch(("tag3", rule2), ("tag4", rule3))
rule3 ::= "c"
"""
    expected = """root ::= TagDispatch(
  ("tag1", rule1),
  ("tag2", rule2),
  ("tag3", rule3),
  loop_after_dispatch=true,
  excludes=()
)
rule1 ::= (("a" rule1_1) | ("zzz"))
rule2 ::= ((rule2_1) | (rule2_2))
rule3 ::= (("c"))
rule1_1 ::= TagDispatch(
  ("tag1", rule2),
  ("tag2", rule3),
  loop_after_dispatch=true,
  excludes=()
)
rule2_1 ::= TagDispatch(
  ("tag1", rule2),
  ("tag2", rule3),
  loop_after_dispatch=true,
  excludes=()
)
rule2_2 ::= TagDispatch(
  ("tag3", rule2),
  ("tag4", rule3),
  loop_after_dispatch=true,
  excludes=()
)
"""
    grammar = xgr.Grammar.from_ebnf(before)
    after = str(grammar)
    assert after == expected


def test_e2e_tag_dispatch_roundtrip():
    """Checks the printed result can be parsed, and the parsing-printing process is idempotent."""
    before = r"""root ::= TagDispatch(
  ("tag1", rule1),
  ("tag2", rule2),
  ("tag3", rule3),
  loop_after_dispatch=false,
  excludes=()
)
rule1 ::= (("a"))
rule2 ::= (("b"))
rule3 ::= (("c"))
"""
    grammar_1 = xgr.Grammar.from_ebnf(before)
    output_string_1 = str(grammar_1)
    grammar_2 = xgr.Grammar.from_ebnf(output_string_1)
    output_string_2 = str(grammar_2)
    assert before == output_string_1
    assert output_string_1 == output_string_2


ebnf_str__expected_error_regex__test_tag_dispatch_parser_errors = [
    (
        'root ::= TagDispatch(("", rule1))\nrule1 ::= "a"',
        "EBNF parser error at line 1, column 21: Tag must be a non-empty string literal",
    ),
    (
        'root ::= TagDispatch(("tag1", undefined_rule))',
        'EBNF parser error at line 1, column 21: Rule "undefined_rule" is not defined',
    ),
    (
        'root ::= TagDispatch("tag1", rule1)',
        "EBNF parser error at line 1, column 21: Each tag dispatch element must be a tuple",
    ),
    (
        'root ::= TagDispatch(("tag1" rule1))',
        "EBNF parser error at line 1, column 30: Expect , or \\) in tuple",
    ),
    (
        'root ::= TagDispatch(("tag1", rule1), stop_str=true)\nrule1 ::= "a"',
        "EBNF parser error at line 1, column 21: Unknown named argument for TagDispatch: stop_str",
    ),
    (
        'root ::= TagDispatch(("tag1", rule1), stop_eos=false)\nrule1 ::= "a"',
        "EBNF parser error at line 1, column 21: Unknown named argument for TagDispatch: stop_eos",
    ),
    (
        'root ::= TagDispatch(("tag1", rule1), excludes=("tag1"))\nrule1 ::= "a"',
        "EBNF parser error at line 1, column 21: Exclude string must not be a prefix of trigger string: tag1",
    ),
]


@pytest.mark.parametrize(
    "ebnf_str, expected_error_regex",
    ebnf_str__expected_error_regex__test_tag_dispatch_parser_errors,
)
def test_tag_dispatch_parser_errors(ebnf_str: str, expected_error_regex: Optional[str]):
    with pytest.raises(RuntimeError, match=expected_error_regex):
        _ebnf_to_grammar_no_normalization(ebnf_str)


def test_trie():
    """Test Trie parsing and printing."""
    before = """root ::= Trie(("abc", "abd", "hello"))
"""
    expected = """root ::= ((Trie(("abc", "abd", "hello"), neg=false)))
"""
    grammar = _ebnf_to_grammar_no_normalization(before)
    assert str(grammar) == expected

    before = """root ::= Trie(("abc", "hello"), neg=true)
"""
    expected = """root ::= ((Trie(("abc", "hello"), neg=true)))
"""
    grammar = _ebnf_to_grammar_no_normalization(before)
    assert str(grammar) == expected


def test_trie_end_to_end():
    grammar = xgr.Grammar.from_ebnf('root ::= Trie(("abc", "abd", "hello"))')
    assert str(grammar) == 'root ::= Trie(("abc", "abd", "hello"), neg=false)\n'
    accepted = ["abc", "abd", "hello"]
    rejected = ["", "ab", "abe", "abcd", "hell", "helloo", "x"]
    for s in accepted:
        assert _is_grammar_accept_string(grammar, s), s
    for s in rejected:
        assert not _is_grammar_accept_string(grammar, s), s


def test_trie_neg_end_to_end():
    # The negated trie rejects any string that has a pattern as a prefix.
    grammar = xgr.Grammar.from_ebnf('root ::= Trie(("abc", "hello"), neg=true)')
    accepted = ["", "ab", "abd", "hell", "x"]
    rejected = ["abc", "hello", "abcd", "helloo"]
    for s in accepted:
        assert _is_grammar_accept_string(grammar, s), s
    for s in rejected:
        assert not _is_grammar_accept_string(grammar, s), s


def test_product():
    """Test Product parsing and printing."""
    before = """root ::= Product(rule1, rule2)
rule1 ::= [a-m] [a-z]
rule2 ::= [k-z] [a-c]
"""
    expected = """root ::= ((Product(rule1, rule2)))
rule1 ::= (([a-m] [a-z]))
rule2 ::= (([k-z] [a-c]))
"""
    grammar = _ebnf_to_grammar_no_normalization(before)
    assert str(grammar) == expected


def test_product_end_to_end():
    # The intersection of rule1 and rule2 is [k-m] [a-c]
    grammar = xgr.Grammar.from_ebnf(
        """root ::= Product(rule1, rule2)
rule1 ::= [a-m] [a-z]
rule2 ::= [k-z] [a-c]
"""
    )
    accepted = ["ka", "mc", "lb"]
    rejected = ["", "ja", "na", "kd", "k", "kab"]
    for s in accepted:
        assert _is_grammar_accept_string(grammar, s), s
    for s in rejected:
        assert not _is_grammar_accept_string(grammar, s), s


def test_product_multiple_rules_with_epsilon():
    # rule1: starts with "a"; rule2: length 2 or 3; rule3: ends with "c"
    grammar = xgr.Grammar.from_ebnf(
        """root ::= Product(rule1, rule2, rule3)
rule1 ::= "a" [a-z]*
rule2 ::= [a-z] [a-z] | [a-z] [a-z] [a-z]
rule3 ::= [a-z]* "c"
"""
    )
    accepted = ["ac", "abc", "axc"]
    rejected = ["", "a", "abcd", "bc", "ab", "abcc"]
    for s in accepted:
        assert _is_grammar_accept_string(grammar, s), s
    for s in rejected:
        assert not _is_grammar_accept_string(grammar, s), s


def test_product_of_tries():
    grammar = xgr.Grammar.from_ebnf(
        """root ::= Product(t1, t2)
t1 ::= Trie(("cat", "dog", "cow"))
t2 ::= Trie(("dog", "cow", "pig"))
"""
    )
    accepted = ["dog", "cow"]
    rejected = ["", "cat", "pig"]
    for s in accepted:
        assert _is_grammar_accept_string(grammar, s), s
    for s in rejected:
        assert not _is_grammar_accept_string(grammar, s), s


def test_product_of_products():
    # A product can reference another product rule; the factor FSMs are built in dependency
    # order and reused directly.
    grammar = xgr.Grammar.from_ebnf(
        """root ::= Product(p1, rule3)
p1 ::= Product(rule1, rule2)
rule1 ::= "a" [a-z]*
rule2 ::= [a-z]* "c"
rule3 ::= [a-z] [a-z]
"""
    )
    # p1: starts with "a" and ends with "c"; rule3: length 2 => only "ac" is accepted
    accepted = ["ac"]
    rejected = ["", "a", "aa", "cc", "abc"]
    for s in accepted:
        assert _is_grammar_accept_string(grammar, s), s
    for s in rejected:
        assert not _is_grammar_accept_string(grammar, s), s


def test_product_allow_empty():
    # The product of two nullable rules accepts the empty string.
    grammar = xgr.Grammar.from_ebnf(
        """root ::= "<" Product(rule1, rule2) ">"
rule1 ::= [ab]*
rule2 ::= [bc]*
"""
    )
    accepted = ["<>", "<b>", "<bb>"]
    rejected = ["<a>", "<c>", "<ab>"]
    for s in accepted:
        assert _is_grammar_accept_string(grammar, s), s
    for s in rejected:
        assert not _is_grammar_accept_string(grammar, s), s


ebnf_str__expected_error_regex__test_trie_product_parser_errors = [
    ('root ::= Trie(("a"), unknown=true)', "Unknown named argument for Trie: unknown"),
    ('root ::= Trie("a", "b")', "Trie\\(\\) requires exactly one tuple of strings"),
    ("root ::= Trie(())", "Trie\\(\\) requires at least one pattern string"),
    ('root ::= Trie((""))', "Trie pattern must be a non-empty string literal"),
    ('root ::= Trie(("a", 1))', "Trie pattern must be a non-empty string literal"),
    ('root ::= Trie(("a"), neg=1)', "neg must be a boolean literal"),
    (
        'root ::= Product(rule1)\nrule1 ::= "a"',
        "Product\\(\\) requires at least two rule references",
    ),
    (
        'root ::= Product(rule1, rule2, x=true)\nrule1 ::= "a"\nrule2 ::= "b"',
        "Product\\(\\) does not accept named arguments",
    ),
    ('root ::= Product(rule1, "b")\nrule1 ::= "a"', "Each argument of Product must be a rule name"),
    ('root ::= Product(rule1, undefined)\nrule1 ::= "a"', 'Rule "undefined" is not defined'),
]


@pytest.mark.parametrize(
    "ebnf_str, expected_error_regex",
    ebnf_str__expected_error_regex__test_trie_product_parser_errors,
)
def test_trie_product_parser_errors(ebnf_str: str, expected_error_regex: Optional[str]):
    with pytest.raises(RuntimeError, match=expected_error_regex):
        _ebnf_to_grammar_no_normalization(ebnf_str)


def test_trie_ac_automaton_product():
    grammar_str = """
    root ::= Product(rule1, rule2) "end"
    rule1 ::= Trie(("cat", "dog", "cow"), neg = true)
    rule2 ::= TagDispatch(excludes = ("end"), loop_after_dispatch = false)
    """
    accepted = ["sheepabcend", "goat12391038end", "pig12i-i9-end"]
    rejected = ["", "catend", "dog12391038end", "cow12i-i9-end", "sheep"]
    grammar = xgr.Grammar.from_ebnf(grammar_str)

    for s in accepted:
        assert _is_grammar_accept_string(grammar, s), s
    for s in rejected:
        assert not _is_grammar_accept_string(grammar, s), s


@pytest.mark.hf_token_required
def test_trie_ac_automaton_product_compilation():
    # 10 patterns of length 10: "abcdefghij", "bcdefghijk", ..., "jklmnopqrs"
    patterns = [string.ascii_lowercase[i : i + 10] for i in range(10)]
    patterns_str = ", ".join(f'"{p}"' for p in patterns)
    grammar_str = f"""root ::= Product(rule1, rule2) "</parameter>"
rule1 ::= Trie(({patterns_str}), neg = true)
rule2 ::= TagDispatch(excludes = ("</parameter>"), loop_after_dispatch = false)
"""
    grammar = xgr.Grammar.from_ebnf(grammar_str)

    model_name = "meta-llama/Llama-3.2-1B-Instruct"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    config = AutoConfig.from_pretrained(model_name)
    tokenizer_info = xgr.TokenizerInfo.from_huggingface(tokenizer, vocab_size=config.vocab_size)
    compiler = xgr.GrammarCompiler(tokenizer_info)

    time_start = time.monotonic_ns()
    compiled_grammar = compiler.compile_grammar(grammar)
    time_end = time.monotonic_ns()
    print(f"Time for compilation: {(time_end - time_start) / 1e6} ms")

    matcher = xgr.GrammarMatcher(compiled_grammar, terminate_without_stop_token=True)
    token_bitmask = xgr.allocate_token_bitmask(1, tokenizer_info.vocab_size)
    matcher.fill_next_token_bitmask(token_bitmask)
    assert matcher.accept_string("some free text</parameter>")
    assert matcher.is_terminated()


if __name__ == "__main__":
    pytest.main(sys.argv)
