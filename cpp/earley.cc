#include "earley.h"

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "fsm.h"
#include "support/logging.h"
#include "support/utils.h"
namespace xgrammar {
Result<RegexIR> RulestoIR(
    const std::vector<std::string>& rules, const std::unordered_map<std::string, int>& rule_map
);
void EarleyParser::Complete(const State& state) {
  auto& state_map = states[state.pos];
  if (state_map.find(state.fsm_num) == state_map.end()) {
    XGRAMMAR_LOG(WARNING) << "State not found in map.";
    return;
  }
  auto& state_set = state_map[state.fsm_num];
  for (const auto& parent : state_set) {
    std::vector<int> dest;
    fsms[parent.fsm_num].fsm.Advance(std::vector<int>{parent.node_num}, state.fsm_num, &dest, true);
    for (const auto& d : dest) {
      process_queue.emplace(d, parent.pos, parent.fsm_num, fsms[parent.fsm_num].IsEndNode(d));
    }
  }
  return;
}

void EarleyParser::PredictAndScan(const State& state, const char& character) {
  auto& fsm = fsms[state.fsm_num];
  std::unordered_set<int> closure;
  fsm.fsm.GetEpsilonClosure(state.node_num, &closure);
  for (const auto& epsilon_state : closure) {
    process_queue.emplace(epsilon_state, state.pos, state.fsm_num, fsm.IsEndNode(epsilon_state));
  }
  std::vector<int> dest;
  fsm.fsm.Advance(std::vector<int>{state.node_num}, character, &dest);
  for (const auto& d : dest) {
    next_states.emplace(d, state.pos + 1, state.fsm_num, fsm.IsEndNode(d));
  }
  std::unordered_set<int> rules;
  fsm.GetPossibleRules(state.node_num, &rules);
  for (const auto& rule : rules) {
    if (states.back().find(rule) == states.back().end()) {
      states.back()[rule] = std::unordered_set<State>();
    }
    states.back()[rule].insert(state);
    process_queue.emplace(
        fsms[rule].StartNode(),
        states.size() - 1,
        rule,
        fsms[rule].IsEndNode(fsms[rule].StartNode())
    );
  }
  return;
}

void EarleyParser::Reset() {
  states.clear();
  current_states.clear();
  next_states.clear();
  process_queue = std::queue<State>();
  current_states.emplace(
      fsms[root_fsm_num].StartNode(),
      0,
      root_fsm_num,
      fsms[root_fsm_num].IsEndNode(fsms[root_fsm_num].StartNode())
  );
  return;
}

void EarleyParser::Consume(const char& character) {
  for (const auto& state : current_states) {
    process_queue.push(state);
  }
  states.emplace_back();
  current_states.clear();
  next_states.clear();
  while (!process_queue.empty()) {
    const auto& state = process_queue.front();
    process_queue.pop();
    if (current_states.find(state) != current_states.end()) {
      continue;
    }
    current_states.insert(state);
    if (state.accept) {
      Complete(state);
    }
    PredictAndScan(state, character);
  }
  current_states = next_states;
  return;
}

bool PythonParser::Parse(const std::string& input) {
  for (const auto& token : input) {
    if (line_start && token == ' ') {
      now_indent++;
      continue;
    }
    if (token == '\n') {
      if (line_start) {
        continue;
      }
      line_start = true;
      last_indent = now_indent;
      now_indent = 0;
      parser.Consume('\n');
      continue;
    }
    if (!line_start) {
      parser.Consume(token);
      continue;
    }
    // line_start is true, and it's the first token which isn't ' ' and '\n'.
    line_start = false;
    if (now_indent % indent_whitespace != 0) {
      XGRAMMAR_LOG(FATAL) << "Indentation error: " << now_indent << " is not a multiple of "
                          << indent_whitespace;
      return false;
    }
    if (now_indent == last_indent) {
      continue;
    }
    if (now_indent > last_indent) {
      int depth = (now_indent - last_indent) / indent_whitespace;
      if (depth > 1) {
        XGRAMMAR_LOG(FATAL) << "Indentation error: " << now_indent << " is not a multiple of "
                            << indent_whitespace;
        return false;
      }
      // TODO: parser.Consume("INDENT");
      continue;
    }
    // now_indent < last_indent
    int depth = (last_indent - now_indent) / indent_whitespace;
    for (int i = 0; i < depth; i++) {
      // TODO: parser.Consume("DEDENT");
    }
  }
  return !parser.current_states.empty();
}

void PythonParser::Reset() {
  parser.Reset();
  line_start = true;
  last_indent = 0;
  now_indent = 0;
  return;
}

Result<PythonParser> GrammarToParser(const std::string& grammar, const int& indent_num) {
  PythonParser parser;
  parser.indent_whitespace = indent_num;
  size_t now_parsing = 0;
  std::vector<std::string> lines;
  std::unordered_map<int, std::vector<std::string>> rules;
  // Build the rule mapping.
  while (grammar.find('\n', now_parsing) != std::string::npos) {
    auto pos = grammar.find('\n', now_parsing);
    lines.push_back(grammar.substr(now_parsing, pos - now_parsing));
    now_parsing = pos + 1;
  }
  if (now_parsing != grammar.size()) {
    lines.push_back(grammar.substr(now_parsing));
  }
  for (const auto& line : lines) {
    if (line.find("::=") == std::string::npos) {
      XGRAMMAR_LOG(FATAL) << "Invalid grammar: " << line;
      return Result<PythonParser>::Err(std::make_shared<Error>("Invalid grammar"));
    }
    if (line.find_first_of("::=") != line.find_last_of("::=")) {
      XGRAMMAR_LOG(FATAL) << "Invalid grammar: " << line;
      return Result<PythonParser>::Err(std::make_shared<Error>("Invalid grammar"));
    }
    auto pos = line.find("::=");
    auto lhs = line.substr(0, pos);
    auto rhs = line.substr(pos + 3);
    while ((!lhs.empty()) && (lhs.back() == ' ')) {
      lhs.pop_back();
    }
    while ((!lhs.empty()) && (lhs.front() == ' ')) {
      lhs.erase(lhs.begin());
    }
    if (parser.rule_map.find(lhs) == parser.rule_map.end()) {
      if (lhs == parser.root_rule) {
        parser.parser.root_fsm_num = parser.rule_map.size();
      }
      parser.rule_map[lhs] = parser.rule_map.size();
      parser.rule_name_map[parser.rule_map[lhs]] = lhs;
      rules[parser.rule_map[lhs]] = std::vector<std::string>({rhs});
    } else {
      rules[parser.rule_map[lhs]].push_back(rhs);
    }
  }
  for (const auto& rule : rules) {
    auto ir = RulestoIR(rule.second, parser.rule_map);
    if (ir.IsErr()) {
      return Result<PythonParser>::Err(ir.UnwrapErr());
    }
    auto built = ir.Unwrap().Build();
    if (built.IsErr()) {
      return Result<PythonParser>::Err(built.UnwrapErr());
    }
    auto fsm = built.Unwrap();
    fsm.SimplifyEpsilon();
    fsm.SimplifyTransition();
    CompactFSMWithStartEnd compact_fsm;
    compact_fsm.fsm = fsm.fsm.ToCompact();
    compact_fsm.ends = std::move(fsm.ends);
    compact_fsm.start = std::move(fsm.start);
    parser.parser.fsms.push_back(std::move(compact_fsm));
  }
  XGRAMMAR_LOG(FATAL) << "Not implemented yet.";
  return Result<PythonParser>::Ok(parser);
};

Result<RegexIR> RulestoIR(
    const std::vector<std::string>& rules, const std::unordered_map<std::string, int>& rule_map
) {
  if (rules.empty()) {
    return Result<RegexIR>::Err(std::make_shared<Error>("Empty rules."));
  }
  if (rules.size() == 1) {
    return RuleToIR(rules[0], rule_map);
  }
  std::vector<RegexIR> irs;
  for (const auto& rule : rules) {
    auto ir = RuleToIR(rule, rule_map);
    if (ir.IsErr()) {
      return ir;
    }
    irs.push_back(ir.Unwrap());
  }
  RegexIR::Union union_node;
  union_node.nodes.reserve(irs.size());
  for (const auto& ir : irs) {
    RegexIR::Bracket bracket;
    for (const auto& node : ir.nodes) {
      bracket.nodes.push_back(node);
    }
    union_node.nodes.push_back(std::move(bracket));
  }
  RegexIR result;
  result.nodes.push_back(std::move(union_node));
  return Result<RegexIR>::Ok(result);
}

Result<RegexIR> RuleToIR(
    const std::string& rule, const std::unordered_map<std::string, int>& rule_map
) {
  // TODO:
  return Result<RegexIR>::Err(std::make_shared<Error>("Not implemented yet."));
}
}  // namespace xgrammar
