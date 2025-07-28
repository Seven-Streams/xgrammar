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
root ::= statements? ENDMARKER

interactive ::= statement_newline

eval ::= expressions NEWLINE* ENDMARKER

func_type ::= "(" type_expressions? [)]  "->" expression NEWLINE* ENDMARKER

statements ::= statement+

statement ::= compound_stmt | simple_stmts

statement_newline ::= (
    compound_stmt NEWLINE
    | simple_stmts
    | NEWLINE
    | ENDMARKER
)
simple_stmts ::= simple_stmt (";" simple_stmt)* ";"? NEWLINE

simple_stmt ::= (
    assignment
    | type_alias
    | star_expressions
    | return_stmt
    | import_stmt
    | raise_stmt
    | "pass"
    | del_stmt
    | yield_stmt
    | assert_stmt
    | "break"
    | "continue"
    | global_stmt
    | nonlocal_stmt
)

compound_stmt ::= (
    function_def
    | if_stmt
    | class_def
    | with_stmt
    | for_stmt
    | try_stmt
    | while_stmt
    | match_stmt
)

assignment ::= (
    NAME ":" expression ("=" annotated_rhs)?
    | ("(" single_target [)]
         | single_subscript_attribute_target) ":" expression ("=" annotated_rhs)?
    | (star_targets "=" )+ (yield_expr | star_expressions) [TYPE_COMMENT]
    | single_target augassign (yield_expr | star_expressions)
)

annotated_rhs ::= yield_expr | star_expressions

augassign ::= (
    "+="
    | "-="
    | "*="
    | "@="
    | "/="
    | "%="
    | "&="
    | "|="
    | "^="
    | "<<="
    | ">>="
    | "**="
    | "//="
)

return_stmt ::= "return" star_expressions?

raise_stmt ::= "raise" (expression ("from" expression )?)?

global_stmt ::= "global" NAME ("," NAME)*

nonlocal_stmt ::= "nonlocal" NAME ("," NAME)*

del_stmt ::= "del" del_targets

yield_stmt ::= yield_expr

assert_stmt ::= "assert" expression ("," expression)?

import_stmt ::= import_name | import_from


import_name ::= "import" dotted_as_names

import_from ::= (
    "from" ("." | "...")* dotted_name "import" import_from_targets
    | "from" ("." | "...")+ "import" import_from_targets
)

import_from_targets ::= (
    "(" import_from_as_names ","? [)]
    | import_from_as_names
    | "*"
)

import_from_as_names ::= import_from_as_name ("," import_from_as_name)*

import_from_as_name ::= NAME ("as" NAME)?

dotted_as_names ::= dotted_as_name ("," dotted_as_name)*

dotted_as_name ::= dotted_name ("as" NAME)

dotted_name ::= (dotted_name ".")? NAME

block ::= (
    NEWLINE INDENT statements DEDENT
    | simple_stmts
)

decorators ::= ("@" named_expression NEWLINE )+

class_def ::= decorators? class_def_raw

class_def_raw ::= "class" NAME (type_params)? ("(" arguments? [)])? ":" block


function_def ::= decorators? function_def_raw

function_def_raw ::= "async"? "def" NAME type_params? "(" params? [)] ("->" expression )? ":" (func_type_comment)? block

params ::= parameters

parameters ::= (
    slash_no_default param_no_default* param_with_default* (star_etc)?
    | slash_with_default param_with_default* (star_etc)?
    | param_no_default+ param_with_default* (star_etc)?
    | param_with_default+ (star_etc)?
    | star_etc
)

slash_no_default ::= param_no_default+ "/" ","?

slash_with_default ::= param_no_default* param_with_default+ "/" ","?

star_etc ::= (
    "*" param_no_default param_maybe_default* kwds?
    | "*" param_no_default_star_annotation param_maybe_default* kwds?
    | "*" "," param_maybe_default+ kwds?
    | kwds
)

kwds ::= "**" param_no_default

param_no_default ::= param ","? TYPE_COMMENT?

param_no_default_star_annotation ::= param_star_annotation ","? TYPE_COMMENT?

param_with_default ::= param default ","? TYPE_COMMENT?

param_maybe_default ::= param default? ","? TYPE_COMMENT?

param ::= NAME annotation?

param_star_annotation ::= NAME star_annotation

annotation ::= ":" expression

star_annotation ::= ":" star_expression

default ::= "=" expression

if_stmt ::= (
    "if" named_expression ":" block elif_stmt
    | "if" named_expression ":" block else_block?
)

elif_stmt ::= (
    "elif" named_expression ":" block elif_stmt
    | "elif" named_expression ":" block else_block?
)

else_block ::= "else" ":" block

while_stmt ::= "while" named_expression ":" block else_block?

for_stmt ::= "async"? "for" star_targets "in" star_expressions ":" TYPE_COMMENT? block else_block?

with_stmt ::= (
    "async"? "with" "(" with_item ("," with_item)*  ","? [)] ":" TYPE_COMMENT? block
    | "async"? "with" with_item ("," with_item)* ":" TYPE_COMMENT? block
)

with_item ::= expression ("as" star_target)?

try_stmt ::= (
    "try" ":" block finally_block
    | "try" ":" block except_block+ (else_block)? (finally_block)?
    | "try" ":" block except_star_block+ (else_block)? (finally_block)?
)

except_block ::= (
    "except" expression ("as" NAME)? ":" block
    | "except" ":" block
)

except_star_block ::= "except" "*" expression ("as" NAME)? ":" block

finally_block ::= "finally" ":" block

match_stmt ::= "match" subject_expr ":" NEWLINE INDENT case_block+ DEDENT

subject_expr ::= (
    star_named_expression "," star_named_expressions?
    | named_expression
)

case_block ::= "case" patterns guard? ":" block

guard ::= "if" named_expression

patterns ::= open_sequence_pattern | pattern

pattern ::= as_pattern | or_pattern

as_pattern ::= or_pattern "as" pattern_capture_target

or_pattern ::= closed_pattern ("|" closed_pattern)*

closed_pattern ::= (
    literal_pattern
    | capture_pattern
    | wildcard_pattern
    | value_pattern
    | group_pattern
    | sequence_pattern
    | mapping_pattern
    | class_pattern
)

literal_pattern ::= (
     signed_number
    | complex_number
    | strings
    | "None"
    | "True"
    | "False"
)

literal_expr ::= (
    signed_number
    | complex_number
    | strings
    | "None"
    | "True"
    | "False"
)

complex_number ::= (
    signed_real_number "+" imaginary_number
    | signed_real_number "-" imaginary_number
)

signed_number ::= "-"? NUMBER

signed_real_number ::= "-"? real_number

real_number ::= NUMBER

imaginary_number ::= NUMBER

capture_pattern ::= pattern_capture_target

pattern_capture_target ::= NAME

wildcard_pattern ::= "_"

value_pattern ::= attr

attr ::= name_or_attr "." NAME

name_or_attr ::= attr | NAME

group_pattern ::= "(" pattern [)]

sequence_pattern ::= "[" maybe_sequence_pattern? "]" | "(" open_sequence_pattern? [)]

open_sequence_pattern ::= maybe_star_pattern "," maybe_sequence_pattern?

maybe_sequence_pattern ::= maybe_star_pattern ("," maybe_star_pattern)* ","?

maybe_star_pattern ::= star_pattern | pattern

star_pattern ::= "*" pattern_capture_target | "*" wildcard_pattern

mapping_pattern ::= (
    "{" "}"
    | "{" double_star_pattern ","? "}"
    | "{" items_pattern "," double_star_pattern ","? "}"
    | "{" items_pattern ","? "}"
)

items_pattern ::= key_value_pattern ("," key_value_pattern)*

key_value_pattern ::= (literal_expr | attr) ":" pattern

double_star_pattern ::= "**" pattern_capture_target

class_pattern ::= (
    name_or_attr "(" [)]
    | name_or_attr "(" positional_patterns ","? [)]
    | name_or_attr "(" keyword_patterns ","? [)]
    | name_or_attr "(" positional_patterns "," keyword_patterns ","? [)]
)

positional_patterns ::= pattern ("," pattern)*

keyword_patterns ::= keyword_pattern ("," keyword_pattern)*

keyword_pattern ::= NAME "=" pattern

type_alias ::= "type" NAME type_params? "=" expression

type_params ::= "[" type_param_seq "]"

type_param_seq ::= type_param ("," type_param)* ","?

type_param ::= (
    NAME type_param_bound? type_param_default?
    | "*" NAME type_param_starred_default?
    | "**" NAME type_param_default?
)

type_param_bound ::= ":" expression

type_param_default ::= "=" expression

type_param_starred_default ::= "=" star_expression

expressions ::= expression ("," expression )* ","?

expression ::= (
    disjunction ("if" disjunction "else" expression)?
    | lambdef
)

yield_expr ::= (
    "yield" "from" expression
    | "yield" star_expressions?
)

star_expressions ::= star_expression ("," star_expression )* ","?

star_expression ::= "*" bitwise_or | expression

star_named_expressions ::= star_named_expression ("," star_named_expression)* ","?

star_named_expression ::= "*" bitwise_or | named_expression

assignment_expression ::= NAME ":=" expression

named_expression ::= assignment_expression | expression

disjunction ::= conjunction ("or" conjunction )*

conjunction ::= inversion ("and" inversion )*

inversion ::= "not" inversion | comparison

comparison ::= bitwise_or compare_op_bitwise_or_pair*

compare_op_bitwise_or_pair ::= (
    eq_bitwise_or
    | noteq_bitwise_or
    | lte_bitwise_or
    | lt_bitwise_or
    | gte_bitwise_or
    | gt_bitwise_or
    | notin_bitwise_or
    | in_bitwise_or
    | isnot_bitwise_or
    | is_bitwise_or
)

eq_bitwise_or ::= "==" bitwise_or

noteq_bitwise_or ::= "!=" bitwise_or

lte_bitwise_or ::= "<=" bitwise_or

lt_bitwise_or ::= "<" bitwise_or

gte_bitwise_or ::= ">=" bitwise_or

gt_bitwise_or ::= ">" bitwise_or

notin_bitwise_or ::= "not" "in" bitwise_or

in_bitwise_or ::= "in" bitwise_or

isnot_bitwise_or ::= "is" "not" bitwise_or

is_bitwise_or ::= "is" bitwise_or

bitwise_or ::= bitwise_or ("|" bitwise_xor)?

bitwise_xor ::= (bitwise_xor "^")? bitwise_and

bitwise_and ::= (bitwise_and "&")? shift_expr

shift_expr ::= (
     shift_expr "<<" sum
    | shift_expr ">>" sum
    | sum
)


sum ::= (
     sum "+" term
    | sum "-" term
    | term
)

term ::= (
    term "*" factor
    | term "/" factor
    | term "//" factor
    | term "%" factor
    | term "@" factor
    | factor
)

factor ::= (
    "+" factor
    | "-" factor
    | "~" factor
    | power
)

power ::= await_primary ("**" factor)?

await_primary ::= "await" primary

primary ::= (
    primary "." NAME
    | primary genexp
    | primary "(" arguments? [)]
    | primary "[" slices "]"
    | atom
)

slices ::= (
    slice
    | (slice | starred_expression) ("," (slice | starred_expression))* ","?
)

slice ::= (
    expression? ":" expression? (":" expression?)?
    | named_expression
)

atom ::= (
    NAME
    | "True"
    | "False"
    | "None"
    | strings
    | NUMBER
    | (tuple | group | genexp)
    | (list | listcomp)
    | (dict | set | dictcomp | setcomp)
    | "..."
)

group ::= "(" (yield_expr | named_expression) [)]

lambdef ::= "lambda" lambda_params? ":" expression

lambda_params ::= lambda_parameters

lambda_parameters ::= (
    lambda_slash_no_default lambda_param_no_default* lambda_param_with_default* lambda_star_etc?
    | lambda_slash_with_default lambda_param_with_default* lambda_star_etc?
    | lambda_param_no_default+ lambda_param_with_default* lambda_star_etc?
    | lambda_param_with_default+ lambda_star_etc?
    | lambda_star_etc
)

lambda_slash_no_default ::= lambda_param_no_default+ "/" ","?

lambda_slash_with_default ::= lambda_param_no_default* lambda_param_with_default+ "/" ","?

lambda_star_etc ::= (
    "*" lambda_param_no_default lambda_param_maybe_default* lambda_kwds?
    | "*" "," lambda_param_maybe_default+ lambda_kwds?
    | lambda_kwds
)

lambda_kwds ::= "**" lambda_param_no_default

lambda_param_no_default ::= lambda_param ","?

lambda_param_with_default ::= lambda_param default ","?

lambda_param_maybe_default ::= lambda_param default? ","?

lambda_param ::= NAME

fstring_middle ::= fstring_replacement_field | FSTRING_MIDDLE

fstring_replacement_field ::= "{" annotated_rhs "="? fstring_conversion? fstring_full_format_spec? "}"

fstring_conversion ::= "!" NAME

fstring_full_format_spec ::= ":" fstring_format_spec*

fstring_format_spec ::= FSTRING_MIDDLE | fstring_replacement_field

fstring ::= FSTRING_START fstring_middle* FSTRING_END

string ::= STRING

strings ::= (fstring|string)+

list ::="[" star_named_expressions? "]"

tuple ::= "(" (star_named_expression "," star_named_expressions?)?[)]

set ::= "{" star_named_expressions "}"

dict ::= "{" [double_starred_kvpairs] "}"

double_starred_kvpairs ::= double_starred_kvpair ("," double_starred_kvpair)* ","?

double_starred_kvpair ::= "**" bitwise_or | kvpair

kvpair ::= expression ":" expression

for_if_clauses ::= for_if_clause+

for_if_clause ::= "async"? "for" star_targets "in" disjunction ("if" disjunction )*

listcomp ::= "[" named_expression for_if_clauses "]"

setcomp ::= "{" named_expression for_if_clauses "}"

genexp ::= "(" ( assignment_expression | expression) for_if_clauses [)]

dictcomp ::= "{" kvpair for_if_clauses "}"

arguments ::= args ","?

args ::= (
    (starred_expression | assignment_expression | expression)
    ("," (starred_expression | assignment_expression | expression))* ("," kwargs)?
    | kwargs
)

kwargs ::= (
    kwarg_or_starred ("," kwarg_or_starred)* "," kwarg_or_double_starred ("," kwarg_or_double_starred)*
    | kwarg_or_starred ("," kwarg_or_starred)*
    | kwarg_or_double_starred ("," kwarg_or_double_starred)*
)

starred_expression ::= "*" expression

kwarg_or_starred ::= NAME "=" expression | starred_expression

kwarg_or_double_starred ::= NAME "=" expression | "**" expression

star_targets ::= star_target ("," star_target )* ","?

star_targets_list_seq ::= star_target ("," star_target)* ","?

star_targets_tuple_seq ::= star_target ("," star_target )* ","?

star_target ::= "*"? target_with_star_atom

target_with_star_atom ::= (
    t_primary "." NAME
    | t_primary "[" slices "]"
    | star_atom
)

star_atom ::= (
    NAME
    | "(" target_with_star_atom [)]
    | "(" star_targets_tuple_seq? [)]
    | "[" star_targets_list_seq? "]"
)

single_target ::= (
    single_subscript_attribute_target
    | NAME
    | "(" single_target [)]
)

single_subscript_attribute_target ::= (
    t_primary "." NAME
    | t_primary "[" slices "]"
)

t_primary ::= (
    t_primary "." NAME
    | t_primary "[" slices "]"
    | t_primary genexp
    | t_primary "(" [arguments] [)]
    | atom
)

t_lookahead ::= "(" | "[" | "."

del_targets ::= del_target ("," del_target)* ","?

del_target ::= (
    t_primary "." NAME
    | t_primary "[" slices "]"
    | del_t_atom
)

del_t_atom ::= (
    NAME
    | "(" del_target [)]
    | "(" del_targets? [)]
    | "[" del_targets? "]"
)

type_expressions ::= (
    expression ("," expression)* "," "*" expression "," "**" expression
    | expression ("," expression)* "," "*" expression
    | expression ("," expression)* "," "**" expression
    | "*" expression "," "**" expression
    | "*" expression
    | "**" expression
    | expression ("," expression)*
)

func_type_comment ::= (
    NEWLINE TYPE_COMMENT
    | TYPE_COMMENT
)

ENDMARKER ::= ""
NEWLINE ::= "\\n"
NAME ::= "[a-zA-Z_][a-zA-Z0-9_]*"
INDENT ::= "INDENT"
DEDENT ::= "DEDENT"
TYPE_COMMENT ::= "#" string
NUMBER ::= integer | floatnumber | imagnumber
FSTRING_MIDDLE ::= "{{" | "}}" | [^\\0{}]
FSTRING_START ::= "f\\""
FSTRING_END ::= ["]
STRING ::= stringliteral | bytesliteral

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
exponent      ::= ("e" | "E") ["+" | "-"] digitpart
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
longbyteschar  ::= [\\x0-[\\]-~])";

}  // namespace xgrammar
