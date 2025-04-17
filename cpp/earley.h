#ifndef XGRAMMAR_EARLEY_H_
#define XGRAMMAR_EARLEY_H_
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "fsm.h"
#include "support/utils.h"
namespace xgrammar {
struct State {
  int node_num;
  int pos;
  int fsm_num;
  bool accept;
  State(int node_num, int pos, int fsm_num, bool accept)
      : node_num(node_num), pos(pos), fsm_num(fsm_num), accept(accept) {}
  State() = default;
  bool operator==(const State& rhs) const {
    return node_num == rhs.node_num && pos == rhs.pos && fsm_num == rhs.fsm_num;
  }
};
}  // namespace xgrammar
namespace std {
template <>
struct hash<xgrammar::State> {
  size_t operator()(const xgrammar::State& state) const {
    return std::hash<int>()(state.node_num) ^ std::hash<int>()(state.pos) ^
           std::hash<int>()(state.fsm_num);
  }
};
}  // namespace std
namespace xgrammar {
class EarleyParser {
 public:
  std::vector<CompactFSMWithStartEnd> fsms;
  std::vector<std::unordered_map<int, std::unordered_set<State>>> states;
  std::unordered_set<State> current_states;
  std::unordered_set<State> next_states;
  std::queue<State> process_queue;
  int root_fsm_num = 0;

  /*!
    \brief EarleyParser Complete function.
    \param state The state which is accepted, and then complete the parent
    states.
  */
  void Complete(const State& state);

  /*!
    \brief EarleyParser Scan and Predict function.
    \param state The state which should be checked.
    \param character The character which should be scanned.
  */
  void PredictAndScan(const State& state, const char& character);

  /*!
    \brief EarleyParser Parse function.
    \param character The character input.
  */
  void Consume(const char& character);

  /*!
    \brief Reset the parser to the initial state.
  */
  void Reset();
  void BuildMap();
};
class PythonParser {
 public:
  EarleyParser parser;
  std::unordered_map<std::string, int> rule_map;
  // Used for debugging.
  std::unordered_map<int, std::string> rule_name_map;
  int indent_whitespace = 4;
  // True if at the beginning of a line.
  // Used for realizing the indentation.
  bool line_start = true;
  // The indentation of the last line.
  int last_indent = 0;
  // The indentation of the current line.
  int now_indent = 0;

  /*!
    \brief Check if the input string is accepted by the parser.
    \param input The input string.
  */
  bool Parse(const std::string& input);

  /*!
    \brief Reset the parser to the initial state.
  */
  void Reset();
};

Result<PythonParser> GrammarToParser(const std::string& grammar, const int& indent_num = 4);
}  // namespace xgrammar

#endif
