#include "earley.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "support/logging.h"
namespace xgrammar {
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

}  // namespace xgrammar
