/*!
 *  Copyright (c) 2024 by Contributors
 * \file xgrammar/grammar.cc
 */

#include <xgrammar/grammar.h>

#include "grammar_functor.h"
#include "grammar_parser.h"
#include "grammar_serializer.h"
#include "json_schema_converter.h"
#include "regex_converter.h"
#include "structural_tag.h"
#include "support/logging.h"

namespace xgrammar {

std::string Grammar::ToString() const { return GrammarPrinter(*this).ToString(); }

Grammar Grammar::FromEBNF(const std::string& ebnf_string, const std::string& root_rule_name) {
  auto grammar = ParseEBNF(ebnf_string, root_rule_name);
  grammar = GrammarNormalizer().Apply(grammar);
  grammar->allow_empty_rule_ids = AllowEmptyRuleAnalyzer::Apply(grammar);
  RepetitionNormalizer::Apply(&grammar);
  GrammarFSMBuilder::Apply(&grammar);
  return grammar;
}

Grammar Grammar::FromJSONSchema(
    const std::string& schema,
    bool any_whitespace,
    std::optional<int> indent,
    std::optional<std::pair<std::string, std::string>> separators,
    bool strict_mode,
    bool print_converted_ebnf
) {
  auto ebnf_string = JSONSchemaToEBNF(schema, any_whitespace, indent, separators, strict_mode);
  if (print_converted_ebnf) {
    XGRAMMAR_LOG(INFO) << "Converted EBNF: " << ebnf_string << std::endl;
  }
  return FromEBNF(ebnf_string);
}

Grammar Grammar::FromRegex(const std::string& regex, bool print_converted_ebnf) {
  auto ebnf_string = RegexToEBNF(regex);
  if (print_converted_ebnf) {
    XGRAMMAR_LOG(INFO) << "Converted EBNF: " << ebnf_string << std::endl;
  }
  return FromEBNF(ebnf_string);
}

Grammar Grammar::FromStructuralTag(
    const std::vector<StructuralTagItem>& tags, const std::vector<std::string>& triggers
) {
  Grammar grammar = StructuralTagToGrammar(tags, triggers);
  grammar->allow_empty_rule_ids = AllowEmptyRuleAnalyzer::Apply(grammar);
  RepetitionNormalizer::Apply(&grammar);
  GrammarFSMBuilder::Apply(&grammar);
  return grammar;
}

// Optimized json grammar for the speed of the grammar matcher
const std::string kJSONGrammarString = R"(
root ::= (
    "{" [ \n\t]* members_and_embrace |
    "[" [ \n\t]* elements_or_embrace
)
value_non_str ::= (
    "{" [ \n\t]* members_and_embrace |
    "[" [ \n\t]* elements_or_embrace |
    "0" fraction exponent |
    [1-9] [0-9]* fraction exponent |
    "-" [0-9] fraction exponent |
    "-" [1-9] [0-9]* fraction exponent |
    "true" |
    "false" |
    "null"
) (= [ \n\t,}\]])
members_and_embrace ::= ("\"" characters_and_colon [ \n\t]* members_suffix | "}") (= [ \n\t,}\]])
members_suffix ::= (
    value_non_str [ \n\t]* member_suffix_suffix |
    "\"" characters_and_embrace |
    "\"" characters_and_comma [ \n\t]* "\"" characters_and_colon [ \n\t]* members_suffix
) (= [ \n\t,}\]])
member_suffix_suffix ::= (
    "}" |
    "," [ \n\t]* "\"" characters_and_colon [ \n\t]* members_suffix
) (= [ \n\t,}\]])
elements_or_embrace ::= (
    "{" [ \n\t]* members_and_embrace elements_rest [ \n\t]* "]" |
    "[" [ \n\t]* elements_or_embrace elements_rest [ \n\t]* "]" |
    "\"" characters_item elements_rest [ \n\t]* "]" |
    "0" fraction exponent elements_rest [ \n\t]* "]" |
    [1-9] [0-9]* fraction exponent elements_rest [ \n\t]* "]" |
    "-" "0" fraction exponent elements_rest [ \n\t]* "]" |
    "-" [1-9] [0-9]* fraction exponent elements_rest [ \n\t]* "]" |
    "true" elements_rest [ \n\t]* "]" |
    "false" elements_rest [ \n\t]* "]" |
    "null" elements_rest [ \n\t]* "]" |
    "]"
)
elements ::= (
    "{" [ \n\t]* members_and_embrace elements_rest |
    "[" [ \n\t]* elements_or_embrace elements_rest |
    "\"" characters_item elements_rest |
    "0" fraction exponent elements_rest |
    [1-9] [0-9]* fraction exponent elements_rest |
    "-" [0-9] fraction exponent elements_rest |
    "-" [1-9] [0-9]* fraction exponent elements_rest |
    "true" elements_rest |
    "false" elements_rest |
    "null" elements_rest
)
elements_rest ::= (
    "" |
    [ \n\t]* "," [ \n\t]* elements
)
characters_and_colon ::= (
    "\"" [ \n\t]* ":" |
    [^"\\\x00-\x1F] characters_and_colon |
    "\\" escape characters_and_colon
) (=[ \n\t]* [\"{[0-9tfn-])
characters_and_comma ::= (
    "\"" [ \n\t]* "," |
    [^"\\\x00-\x1F] characters_and_comma |
    "\\" escape characters_and_comma
) (=[ \n\t]* "\"")
characters_and_embrace ::= (
    "\"" [ \n\t]* "}" |
    [^"\\\x00-\x1F] characters_and_embrace |
    "\\" escape characters_and_embrace
) (=[ \n\t]* [},])
characters_item ::= (
    "\"" |
    [^"\\\x00-\x1F] characters_item |
    "\\" escape characters_item
) (= [ \n\t]* [,\]])
escape ::= ["\\/bfnrt] | "u" [A-Fa-f0-9] [A-Fa-f0-9] [A-Fa-f0-9] [A-Fa-f0-9]
fraction ::= "" | "." [0-9] [0-9]*
exponent ::= "" |  "e" sign [0-9] [0-9]* | "E" sign [0-9] [0-9]*
sign ::= "" | "+" | "-"
)";

Grammar Grammar::BuiltinJSONGrammar() {
  static const Grammar grammar = FromEBNF(kJSONGrammarString);
  return grammar;
}

Grammar Grammar::Union(const std::vector<Grammar>& grammars) {
  return GrammarUnionFunctor::Apply(grammars);
}

Grammar Grammar::Concat(const std::vector<Grammar>& grammars) {
  return GrammarConcatFunctor::Apply(grammars);
}

std::ostream& operator<<(std::ostream& os, const Grammar& grammar) {
  os << grammar.ToString();
  return os;
}

const std::string KPythonGrammarString = R"(
root ::= statements? ENDMARKER

WB ::= [ ]*

MUSTWB ::= [ ]+

interactive ::= statement_newline

eval ::= expressions NEWLINE* ENDMARKER

func_type ::= "(" WB type_expressions? [)] WB "->" WB expression NEWLINE* ENDMARKER

statements ::= statement+

statement ::= compound_stmt | simple_stmts

statement_newline ::= (
    compound_stmt NEWLINE |
    simple_stmts |
    NEWLINE |
    ENDMARKER
)
simple_stmts ::= simple_stmt (";" WB simple_stmt)* (";" WB)? NEWLINE

simple_stmt ::= (
    assignment |
    type_alias |
    star_expressions |
    return_stmt |
    import_stmt |
    raise_stmt |
    "pass" |
    del_stmt |
    yield_stmt |
    assert_stmt |
    "break" |
    "continue" |
    global_stmt |
    nonlocal_stmt
)

compound_stmt ::= (
    function_def |
    if_stmt |
    class_def |
    with_stmt |
    for_stmt |
    try_stmt |
    while_stmt |
    match_stmt
)

assignment ::= (
    NAME ":" WB expression ("=" WB (yield_expr | star_expressions))? |
    ("(" WB single_target [)] WB |
    single_subscript_attribute_target) ":" WB expression ("=" WB (yield_expr | star_expressions))? |
    (star_targets "=" WB)+ (yield_expr | star_expressions) TYPE_COMMENT? |
    single_target (( [+\-*/@%&|^] | ">>" | "<<" | "//" | "**" ) "=" WB) (yield_expr | star_expressions)
)

return_stmt ::= "return" MUSTWB star_expressions?

raise_stmt ::= "raise" MUSTWB (expression (MUSTWB "from" MUSTWB expression )?)?

global_stmt ::= "global" MUSTWB NAME ("," WB NAME)*

nonlocal_stmt ::= "nonlocal" MUSTWB NAME ("," WB NAME)*

del_stmt ::= "del" MUSTWB del_targets

yield_stmt ::= yield_expr

assert_stmt ::= "assert" MUSTWB expression ("," WB expression)?

import_stmt ::= (
    "import" MUSTWB dotted_as_names |
    "from" MUSTWB ("." | "...")* dotted_name MUSTWB "import" MUSTWB import_from_targets |
    "from" MUSTWB ("." | "...")+ MUSTWB "import" MUSTWB import_from_targets
)

import_from_targets ::= (
    "(" WB import_from_as_names ","? WB [)] WB |
    import_from_as_names |
    "*" WB
)

import_from_as_names ::= NAME (MUSTWB "as" MUSTWB NAME)? ("," WB NAME (MUSTWB "as" MUSTWB NAME)? )*

dotted_as_names ::= dotted_as_name ("," WB dotted_as_name)*

dotted_as_name ::= dotted_name (MUSTWB "as" MUSTWB NAME)

dotted_name ::= (dotted_name "." WB)? NAME

block ::= (
    NEWLINE INDENT statements DEDENT |
    simple_stmts
)

decorators ::= ("@" WB named_expression NEWLINE )+

class_def ::= decorators? "class" MUSTWB NAME (type_params)? ("(" WB arguments? [)] WB)? ":" WB block

function_def ::= decorators? ("async" MUSTWB)? "def" MUSTWB NAME type_params? "(" WB parameters? [)] WB ("->" WB expression )? ":" WB (func_type_comment)? block

parameters ::= (
    param_no_default+ "/" WB ","? WB param_no_default* param_with_default* (star_etc)? |
    param_no_default* param_with_default+ "/" WB ","? WB param_with_default* (star_etc)? |
    param_no_default+ param_with_default* (star_etc)? |
    param_with_default+ (star_etc)? |
    star_etc
)

star_etc ::= (
    "*" WB param_no_default param_maybe_default* kwds? |
    "*" WB NAME ":" WB star_expression ","? WB TYPE_COMMENT? param_maybe_default* kwds? |
    "*" WB "," WB param_maybe_default+ kwds? |
    kwds
)

kwds ::= "**" WB param_no_default

param_no_default ::= NAME WB expression? ","? WB TYPE_COMMENT?

param_with_default ::= NAME WB expression? ","? WB TYPE_COMMENT?

param_maybe_default ::= NAME WB expression? ","? WB TYPE_COMMENT?

default ::= "=" WB expression

if_stmt ::= (
    "if" MUSTWB named_expression ":" WB block elif_stmt |
    "if" MUSTWB named_expression ":" WB block else_block?
)

elif_stmt ::= (
    "elif" MUSTWB named_expression ":" WB block elif_stmt |
    "elif" MUSTWB named_expression ":" WB block else_block?
)

else_block ::= "else" WB ":" WB block

while_stmt ::= "while" MUSTWB named_expression ":" WB block else_block?

for_stmt ::= ("async" MUSTWB)? "for" MUSTWB star_targets MUSTWB "in" MUSTWB star_expressions WB ":" WB TYPE_COMMENT? block else_block?

with_stmt ::= (
    ("async" MUSTWB)? "with" MUSTWB "(" WB with_item ("," WB with_item)*  ","? WB [)] WB ":" WB TYPE_COMMENT? block |
    ("async" MUSTWB)? "with" MUSTWB with_item ("," WB with_item)* ":" WB TYPE_COMMENT? block
)

with_item ::= expression (MUSTWB "as" MUSTWB star_target)?

try_stmt ::= (
    "try" WB ":" WB block finally_block |
    "try" WB ":" WB block except_block+ (else_block)? (finally_block)? |
    "try" WB ":" WB block except_star_block+ (else_block)? (finally_block)?
)

except_block ::= "except" MUSTWB (expression (MUSTWB "as" MUSTWB NAME)?)? ":" WB block

except_star_block ::= "except" MUSTWB "*" WB expression (MUSTWB "as" MUSTWB NAME)? ":" WB block

finally_block ::= "finally" WB ":" WB block

match_stmt ::= "match" MUSTWB (
    star_named_expression "," WB star_named_expressions? |
    named_expression)
 ":" WB NEWLINE INDENT case_block+ DEDENT

case_block ::= "case" MUSTWB patterns ("if" MUSTWB named_expression)? ":" WB block

patterns ::= open_sequence_pattern | pattern

pattern ::= as_pattern | or_pattern

as_pattern ::= or_pattern MUSTWB "as" MUSTWB pattern_capture_target

or_pattern ::= closed_pattern ("|" WB closed_pattern)*

closed_pattern ::= (
    literal_pattern |
    capture_pattern |
    wildcard_pattern |
    value_pattern |
    group_pattern |
    sequence_pattern |
    mapping_pattern |
    class_pattern
)

literal_pattern ::= (
    signed_number |
    complex_number |
    strings |
    "None" WB|
    "True" WB|
    "False" WB
)

complex_number ::= (
    signed_number "+" WB NUMBER |
    signed_number "-" WB NUMBER
)

signed_number ::= ("-" WB)? NUMBER

capture_pattern ::= pattern_capture_target

pattern_capture_target ::= NAME

wildcard_pattern ::= "_" WB

value_pattern ::= attr

attr ::= name_or_attr "." WB NAME

name_or_attr ::= attr | NAME

group_pattern ::= "(" WB pattern [)] WB

sequence_pattern ::= "[" WB maybe_sequence_pattern? "]" WB | "(" WB open_sequence_pattern? [)] WB

open_sequence_pattern ::= maybe_star_pattern "," WB maybe_sequence_pattern?

maybe_sequence_pattern ::= maybe_star_pattern ("," WB maybe_star_pattern)* ","? WB

maybe_star_pattern ::= star_pattern | pattern

star_pattern ::= "*" WB pattern_capture_target | "*" WB wildcard_pattern

mapping_pattern ::= (
    "{" WB "}" WB |
    "{" WB double_star_pattern ","? WB "}" WB |
    "{" WB items_pattern "," double_star_pattern ","? WB "}" WB|
    "{" WB items_pattern ","? WB "}" WB
)

items_pattern ::= key_value_pattern ("," WB key_value_pattern)*

key_value_pattern ::= (literal_pattern | attr) ":" WB pattern

double_star_pattern ::= "**" WB pattern_capture_target

class_pattern ::= (
    name_or_attr "(" WB [)] WB|
    name_or_attr "(" WB positional_patterns ","? WB [)] WB|
    name_or_attr "(" WB keyword_patterns ","? WB [)] WB |
    name_or_attr "(" WB positional_patterns "," WB keyword_patterns ","? WB [)] WB
)

positional_patterns ::= pattern ("," WB pattern)*

keyword_patterns ::= keyword_pattern ("," WB keyword_pattern)*

keyword_pattern ::= NAME "=" WB pattern

type_alias ::= "type" MUSTWB NAME type_params? "=" WB expression

type_params ::= "[" WB type_param_seq "]" WB

type_param_seq ::= type_param ("," WB type_param)* ","? WB

type_param ::= (
    NAME type_param_bound? type_param_default? |
    "*" WB NAME type_param_starred_default? |
    "**" WB NAME type_param_default?
)

type_param_bound ::= ":" WB expression

type_param_default ::= "=" WB expression

type_param_starred_default ::= "=" WB star_expression

expressions ::= expression ("," WB expression )* ","? WB

expression ::= (
    disjunction ("if" MUSTWB disjunction "else" MUSTWB expression)? |
    lambdef
)

yield_expr ::= (
    "yield" MUSTWB "from" MUSTWB expression |
    "yield" MUSTWB star_expressions?
)

star_expressions ::= star_expression ("," WB star_expression )* ","? WB

star_expression ::= "*" WB bitwise_or | expression

star_named_expressions ::= star_named_expression ("," WB star_named_expression)* ","? WB

star_named_expression ::= "*" WB bitwise_or | named_expression

assignment_expression ::= NAME ":=" WB expression

named_expression ::= assignment_expression | expression

disjunction ::= conjunction (MUSTWB "or" MUSTWB conjunction )*

conjunction ::= inversion (MUSTWB "and" MUSTWB inversion )*

inversion ::= "not" MUSTWB inversion | comparison

comparison ::= bitwise_or compare_op_bitwise_or_pair*

compare_op_bitwise_or_pair ::= (
    [!>=<] "=" | [><] | ("not" MUSTWB)? "in" MUSTWB | "is" MUSTWB ("not" MUSTWB)?
) bitwise_or

bitwise_or ::= bitwise_xor ("|" WB bitwise_xor)?

bitwise_xor ::= (bitwise_xor "^" WB)? bitwise_and

bitwise_and ::= (bitwise_and "&" WB)? shift_expr

shift_expr ::= (
     shift_expr "<<" WB sum |
     shift_expr ">>" WB sum |
     sum
)

sum ::= (
     sum "+" WB term |
     sum "-" WB term |
     term
)

term ::= (
    term "*" WB factor |
    term "/" ("/")? WB factor |
    term "%" WB factor |
    term "@" WB factor |
    factor
)

factor ::= (
    "+" WB factor |
    "-" WB factor |
    "~" WB factor |
    power
)

power ::= await_primary ("**" WB factor)?

await_primary ::= ("await" MUSTWB)? primary

primary ::= (
    primary "." WB NAME |
    primary genexp |
    primary "(" WB arguments? [)] WB |
    primary "[" WB slices "]" WB |
    atom
)

slices ::= (
    slice |
    (slice | starred_expression) ("," WB (slice | starred_expression))* ","? WB
)

slice ::= (
    expression? ":" WB expression? (":" WB expression?)? |
    named_expression
)

atom ::= (
    NAME |
    "True" WB|
    "False" WB|
    "None" WB|
    strings |
    NUMBER |
    (tuple | group | genexp) |
    (list | listcomp) |
    (dict | set | dictcomp | setcomp) |
    "..."
)

group ::= "(" WB (yield_expr | named_expression) [)] WB

lambdef ::= "lambda" MUSTWB lambda_params? ":" WB expression

lambda_params ::= lambda_parameters

lambda_parameters ::= (
    lambda_slash_no_default lambda_param_no_default* lambda_param_with_default* lambda_star_etc? |
    lambda_slash_with_default lambda_param_with_default* lambda_star_etc? |
    lambda_param_no_default+ lambda_param_with_default* lambda_star_etc? |
    lambda_param_with_default+ lambda_star_etc? |
    lambda_star_etc
)

lambda_slash_no_default ::= lambda_param_no_default+ "/" WB ","? WB

lambda_slash_with_default ::= lambda_param_no_default* lambda_param_with_default+ "/" WB ","? WB

lambda_star_etc ::= (
    "*" WB lambda_param_no_default lambda_param_maybe_default* lambda_kwds? |
    "*" WB "," WB lambda_param_maybe_default+ lambda_kwds? |
    lambda_kwds
)

lambda_kwds ::= "**" WB lambda_param_no_default

lambda_param_no_default ::= lambda_param ","? WB

lambda_param_with_default ::= lambda_param default ","? WB

lambda_param_maybe_default ::= lambda_param default? ","? WB

lambda_param ::= NAME

fstring_middle ::= fstring_replacement_field | FSTRING_MIDDLE

fstring_replacement_field ::= "{" WB (yield_expr | star_expressions) "="? WB fstring_conversion? fstring_full_format_spec? "}" WB

fstring_conversion ::= "!" WB NAME

fstring_full_format_spec ::= ":" WB fstring_format_spec*

fstring_format_spec ::= FSTRING_MIDDLE | fstring_replacement_field

fstring ::= FSTRING_START fstring_middle* FSTRING_END

string ::= STRING

strings ::= (fstring|string)+

list ::="[" WB star_named_expressions? "]" WB

tuple ::= "(" WB (star_named_expression "," WB star_named_expressions?)? [)] WB

set ::= "{" WB star_named_expressions "}" WB

dict ::= "{" WB double_starred_kvpairs? "}" WB

double_starred_kvpairs ::= double_starred_kvpair ("," WB double_starred_kvpair)* ","? WB

double_starred_kvpair ::= "**" WB bitwise_or | kvpair

kvpair ::= expression ":" WB expression

for_if_clauses ::= for_if_clause+

for_if_clause ::= ("async" MUSTWB)? "for" MUSTWB star_targets MUSTWB "in" MUSTWB disjunction (MUSTWB "if" MUSTWB disjunction )*

listcomp ::= "[" WB named_expression for_if_clauses "]" WB

setcomp ::= "{" WB named_expression for_if_clauses "}" WB

genexp ::= "(" WB (assignment_expression | expression) for_if_clauses [)] WB

dictcomp ::= "{" WB kvpair for_if_clauses "}" WB

arguments ::= args ","? WB

args ::= (
    (starred_expression | assignment_expression | expression)
    ("," WB (starred_expression | assignment_expression | expression))* ("," WB kwargs)? |
    kwargs
)

kwargs ::= (
    kwarg_or_starred ("," WB kwarg_or_starred)* "," WB kwarg_or_double_starred ("," WB kwarg_or_double_starred)* |
    kwarg_or_starred ("," WB kwarg_or_starred)* |
    kwarg_or_double_starred ("," WB kwarg_or_double_starred)*
)

starred_expression ::= "*" WB expression

kwarg_or_starred ::= NAME "=" WB expression | starred_expression

kwarg_or_double_starred ::= NAME "=" WB expression | "**" WB expression

star_targets ::= star_target ("," WB star_target )* ","? WB

star_targets_list_seq ::= star_target ("," WB star_target)* ","? WB

star_targets_tuple_seq ::= star_target ("," WB star_target )* ","? WB

star_target ::= ("*" WB)? target_with_star_atom

target_with_star_atom ::= (
    t_primary "." WB NAME |
    t_primary "[" WB slices "]" WB |
    star_atom
)

star_atom ::= (
    NAME |
    "(" WB target_with_star_atom [)] WB |
    "(" WB star_targets_tuple_seq? [)] WB|
    "[" WB star_targets_list_seq? "]" WB
)

single_target ::= (
    single_subscript_attribute_target |
    NAME |
    "(" WB single_target [)] WB
)

single_subscript_attribute_target ::= (
    t_primary "." WB NAME |
    t_primary "[" WB slices "]" WB
)

t_primary ::= (
    t_primary "." WB NAME |
    t_primary "[" WB slices "]" WB |
    t_primary genexp |
    t_primary "(" WB arguments? [)] WB |
    atom
)

t_lookahead ::= ("(" | "[" | ".") WB

del_targets ::= del_target ("," WB del_target)* ","? WB

del_target ::= (
    t_primary "." WB NAME |
    t_primary "[" WB slices "]" WB |
    del_t_atom
)

del_t_atom ::= (
    NAME |
    "(" WB del_target [)] WB |
    "(" WB del_targets? [)] WB |
    "[" WB del_targets? "]" WB
)

type_expressions ::= (
    expression ("," WB expression)* "," WB "*" WB expression ("," WB "**" WB expression)? |
    expression ("," WB expression)* "," WB "**" WB expression |
    "*" WB expression "," WB "**" WB expression |
    "*" WB expression |
    "**" WB expression |
    expression ("," WB expression)*
)

func_type_comment ::= (
    NEWLINE TYPE_COMMENT |
    TYPE_COMMENT
)

ENDMARKER ::= ""
NEWLINE ::= "\\n"
NAME ::= [a-zA-Z_][a-zA-Z0-9_]* WB
INDENT ::= "INDENT"
DEDENT ::= "DEDENT"
TYPE_COMMENT ::= "#" longstringitem
NUMBER ::= (integer | floatnumber | imagnumber) WB
FSTRING_MIDDLE ::= "{{" | "}}" | [^\\0{}]
FSTRING_START ::= "f\\""
FSTRING_END ::= ["]
STRING ::= (stringliteral | bytesliteral) WB

integer      ::= decinteger | bininteger | octinteger | hexinteger
decinteger   ::= nonzerodigit ("_"? digit)* | "0"+ ("_"? "0")*
bininteger   ::= "0" [bB] ("_"? bindigit)+
octinteger   ::= "0" [oO] ("_"? octdigit)+
hexinteger   ::= "0" [xX] ("_"? hexdigit)+
nonzerodigit ::= [1-9]
digit        ::= [0-9]
bindigit     ::= [01]
octdigit     ::= [0-7]
hexdigit     ::= [0-9a-fA-F]
floatnumber   ::= pointfloat | exponentfloat
pointfloat    ::= [digitpart] fraction | digitpart "."
exponentfloat ::= (digitpart | pointfloat) exponent
digitpart     ::= digit ("_"? digit)*
fraction      ::= "." digitpart
exponent      ::= ("e" | "E") ("+" | "-")? digitpart
imagnumber ::= (floatnumber | digitpart) ("j" | "J")
stringliteral   ::= stringprefix? (shortstring | longstring)
stringprefix    ::= "r" | "u" | "R" | "U" | "f" | "F"
                    | "fr" | "Fr" | "fR" | "FR" | "rf" | "rF" | "Rf" | "RF"
shortstring     ::= "'" shortstringitem* "'" | ["] shortstringitem* ["]
longstring      ::= "'''" longstringitem* "'''" | ["]["]["] longstringitem* ["]["]["]
shortstringitem ::= shortstringchar | stringescapeseq
longstringitem  ::= longstringchar | stringescapeseq
bytesliteral   ::= bytesprefix(shortbytes | longbytes)
bytesprefix    ::= "b" | "B" | "br" | "Br" | "bR" | "BR" | "rb" | "rB" | "Rb" | "RB"
shortbytes     ::= "'" shortbytesitem* "'" | ["] shortbytesitem* ["]
longbytes      ::= "'''" longbytesitem* "'''" | ["]["]["] longbytesitem* ["]["]["]
shortbytesitem ::= shortbyteschar | bytesescapeseq
longbytesitem  ::= longbyteschar | bytesescapeseq
shortstringchar ::= [^\"\\n]
stringescapeseq ::= "\" [\\x0-\\x7f]
longstringchar  ::= [^\\\\]
shortbyteschar ::= [\\x0-\\x09\\x11-!#-[\\]-~]
bytesescapeseq ::= [\\\\] [\\x0-\\x7f]
longbyteschar  ::= [\\x0-[\\]-~]
)";

}  // namespace xgrammar
