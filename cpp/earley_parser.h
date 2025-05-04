/*!
 *  Copyright (c) 2025 by Contributors
 * \file xgrammar/earley_parser.h
 * \brief The header for the definition of the Earley parser.
 */
#ifndef XGRAMMAR_EARLEY_PARSER_H_
#define XGRAMMAR_EARLEY_PARSER_H_
#include <cstdint>
#include <map>
#include <ostream>
#include <queue>
#include <utility>
#include <vector>

#include "grammar_data_structure.h"
#include "support/csr_array.h"
#include "support/utils.h"
#include "xgrammar/grammar.h"

namespace xgrammar {

/* \brief The state of the Earley parser.
  In the implementation, a rule can only be a kchoices or a ktagdispatch.
  A kchoices rule must be composed of some ksequence rules, or a kemptyrule.
  In the ksequence, every element in the sequence must be a kbytestring, a
  kcharacterclass, a kcharacterclassstar, or a rule reference.

  -rule_id: The id of the rule.
  -sequence_id: The id of the sequence in the rule.
  -element_id: The id of the element in the sequence, or the id of the node in
  the tag dispatch fsm.
  -parent_pos: The id of the parent node in the Earley parser. i.e. the rule
  is predicted from the k-th character.
  -sub_element_id: The id of the sub element in the current element, i.e.:
    - kbytestring: the id of the byte in the string.
    - kcharacterclass: How many bytes are left to be read in the utf8 character.
    - kcharacterclassstar: How many bytes are left to be read in the utf8 character.
*/
struct ParserState {
  /*! \brief A sequence_id value of kUnexpandedRuleStartSequenceId means a rule hasn't been
   * expanded.*/
  static constexpr int32_t kUnexpandedRuleStartSequenceId = 128000;

  /*! \brief A parent_id value of kNoParent means this ParserState is the root of the parsing stack.
   */
  static constexpr int32_t kNoPrevInputPos = -1;

  /*! \brief A sequence_id value of kInvalid means the ParserState is invalid.*/
  static constexpr int32_t kInvalid = -2;

  /*! \brief The rule's id.*/
  int32_t rule_id = -1;

  /*! \brief Which choice in this rule is selected. */
  int32_t sequence_id = -1;

  /*! \brief Which element of the choice sequence is to be visited. When the current sequence is
   * a tag dispatch rule, this element id the currently visited node. */
  int32_t element_id = -1;

  /*! \brief The position of the state, i.e. from which position, the rule starts.*/
  int32_t input_pos = -1;

  /*! \brief The id of the sub element in the current selement of the sequence. */
  int32_t sub_element_id = 0;

  /*! \brief If the rule is completed.*/
  bool completed = false;

  constexpr ParserState() = default;

  constexpr ParserState(const ParserState&) = default;

  ParserState& operator=(const ParserState&) = default;

  constexpr ParserState(
      const int32_t& rule_id,
      const int32_t& sequence_id,
      const int32_t& element_id,
      const int32_t& parent_pos,
      const int32_t& sub_element_id
  )
      : rule_id(rule_id),
        sequence_id(sequence_id),
        element_id(element_id),
        input_pos(parent_pos),
        sub_element_id(sub_element_id) {}

  // The element is invalid when sequence_id is -1.
  bool IsInvalid() const { return sequence_id == -1; }

  static ParserState GetInvalidState() { return {-1, -1, -1, -1, -1}; }

  bool operator==(const ParserState& other) const {
    return rule_id == other.rule_id && sequence_id == other.sequence_id &&
           element_id == other.element_id && sub_element_id == other.sub_element_id;
  }

  friend std::ostream& operator<<(std::ostream& os, const ParserState& state) {
    os << "ParserState(rule_id=" << state.rule_id << ", sequence_id=" << state.sequence_id
       << ", element_id=" << state.element_id << ", parent_pos=" << state.input_pos
       << ", sub_element_id=" << state.sub_element_id << ")";
    return os;
  }
};

/*
  When getting the mask of the state, we don't need to consider the parent_pos.
*/
class StateHash {
 public:
  size_t operator()(const ParserState& state) const {
    return HashCombine(state.rule_id, state.sequence_id, state.element_id, state.sub_element_id);
  }
};

/*
  When matching the state, we need to consider the parent_pos, since if two states
  don't have the same parent_pos, they are not the same state.
*/
class StateEqual {
 public:
  bool operator()(const ParserState& lhs, const ParserState& rhs) const {
    return lhs.rule_id == rhs.rule_id && lhs.sequence_id == rhs.sequence_id &&
           lhs.element_id == rhs.element_id && lhs.input_pos == rhs.input_pos &&
           lhs.sub_element_id == rhs.sub_element_id;
  }
};

class StateHashChecker {
 public:
  size_t operator()(const ParserState& state) const {
    return HashCombine(
        state.rule_id, state.sequence_id, state.element_id, state.input_pos, state.sub_element_id
    );
  }
};
class EarleyParser {
 protected:
  /*! \brief The grammar to be parsed. */
  Grammar grammar_;

  /*! \brief rule_id_to_completeable_states[i][j] is the i pos j rule_id states. It's used for
   * completion. */
  std::vector<std::multimap<int32_t, ParserState>> rule_id_to_completeable_states;

  /*!
      \brief The states history. state_stack[i] is a vector storing the states after accepting the
     input[i-1].
   */
  CSRArray<ParserState> scanable_state_history_;

  /*! \brief A temperate vector only used in Advance, used to add states in the
   * scanable_state_history. */
  std::vector<ParserState> tmp_states_to_be_added_;

  /*! \brief It's the processing queue of the earley parser.*/
  std::queue<ParserState> tmp_process_state_queue_;

  /*! The vector to check if a state has been added into the queue.*/
  std::vector<ParserState> tmp_states_visited_in_queue_;

  /*!
  \brief Check if the state has been added into the queue.
  \param state The state to check.
  \return True if in the vector, false otherwise.
*/
  bool IsStateVisitedInQueue(const ParserState& state) const {
    return (
        std::find_if(
            tmp_states_visited_in_queue_.begin(),
            tmp_states_visited_in_queue_.end(),
            [&](const ParserState& s) { return StateEqual()(state, s); }
        ) != tmp_states_visited_in_queue_.end()
    );
  }

  /*!
    \brief Push the state into the queue. If the state is already in the queue,
    then we don't need to push it again.
    \param state The state to be pushed.
  */
  void Enqueue(const ParserState& state) {
    if (!IsStateVisitedInQueue(state)) {
      tmp_process_state_queue_.push(state);
      tmp_states_visited_in_queue_.push_back(state);
    }
    return;
  }

  /*!
    \brief The scanning operation of the Earley parser.
  */
  void Scan(const ParserState& state, const uint8_t& ch);

  /*!
      \brief The completion operation of the Earley parser.
      \return If the state is scannable then we return false(which means that
      the state can't predict and complete), otherwise we return true.
      \details The reason is that if the state can't be scanned, then
      add it into the next states is useless. Moreover, the end
      of the grammar is used to check if the grammar is completed,
      so it should be added into the next states.
  */
  bool Complete(const ParserState& state);

  /*!
      \brief The prediction operation of the Earley parser.
      \return Fitst: If the state scanable, or the state is the end of the grammar,
      then return true, otherwise return false.
      \return Second: If the state is obviously that it can't complete,
      then return false, true otherwise.
  */
  std::pair<bool, bool> Predict(const ParserState& state);

  /*!
    \brief Check if the state is the end of the grammar.
    \param state The state to be checked.
    \return True if the state is the end of the grammar, false otherwise.
  */
  bool IsStateCompleted(const ParserState& state) const;

  /*!
    \brief Check if a character can be accepted.
  */
  bool IsAccepted(const ParserState& state, uint8_t ch) const;

  /*!
    \brief Handle the unexpanded rule, used for pushing initial state.
    \param state The state to be handled.
    \return True if the rule is unexpanded, false otherwise.
  */
  bool ExpandAndEnqueueUnexpandedState(const ParserState& state);

  /*!
    \brief Expand the rule, used for RuleRef and kTagDispatch.
    \param state The state to be expanded, which is the parent state.
    The type of the state is kTagDispatch or kSequence. Moreover, the
    element of the sequence should be a rule reference; the node in
    the kTagDispatch should be an end node.
  */
  void ExpandNextRuleRefElement(const ParserState& state, const Grammar::Impl::RuleExpr& rule_expr);

  /*!
    \brief Advance the parser to the next state, with the sub rule is kCharacterClass.
    \param state The state to be advanced.
    \param ch The character to be advanced.
    \param cur_sequence The sequence of the current state.
    \param sub_sequence The sub sequence to be checked.
    \return The next state, Invalid state if the character is not accepted.
  */
  void AdvanceCharacterClass(
      const ParserState& state,
      const uint8_t& ch,
      const Grammar::Impl::RuleExpr& cur_sequence,
      const Grammar::Impl::RuleExpr& sub_sequence
  );

  /*!
    \brief Advance the parser to the next state, with the sub rule is kByteString.
    \param state The state to be advanced.
    \param ch The character to be advanced.
    \param cur_sequence The sequence of the current state.
    \param sub_sequence The sub sequence to be checked.
    \return The next state, Invalid state if the character is not accepted.
  */
  void AdvanceByteString(
      const ParserState& state,
      const uint8_t& ch,
      const Grammar::Impl::RuleExpr& cur_sequence,
      const Grammar::Impl::RuleExpr& sub_sequence
  );

  /*!
    \brief Advance the parser to the next state, with the sub rule is kCharacterClassStar.
    \param state The state to be advanced.
    \param ch The character to be advanced.
    \param cur_sequence The sequence of the current state.
    \param sub_sequence The sub sequence to be checked.
    \return The next state, Invalid state if the character is not accepted.
  */
  void AdvanceCharacterClassStar(
      const ParserState& state,
      const uint8_t& ch,
      const Grammar::Impl::RuleExpr& cur_sequence,
      const Grammar::Impl::RuleExpr& sub_sequence
  );

 public:
  /*!
   \brief Constructor of the Earley parser.
   \param grammar The grammar to be parsed.
   \param initial_state The initial state to be pushed into the parser.
 */
  EarleyParser(
      const Grammar& grammar, const ParserState& initial_state, const bool& need_expand = true
  );

  /*!
    \brief From the current states, advance to the next state.
    \param ch The character to be advanced.
    \return True if the character is accepted, false otherwise.
    \note If the character isn't accepted, then the states won't be changed.
  */
  bool Advance(const uint8_t& ch);

  /*!
    \brief Remove the newly added states.
    \param count The number of states to be removed.
  */
  void PopLastStates(int32_t count = 1);

  /*!
    \brief Check whether any of the multiple states stored in the parser has already completed.
    \note Since the parser contains multiple parallel states, some may have already completed, while
    others might still be able to accept more characters. \return True if the root rule is
    completed, false otherwise.
  */
  bool IsCompleted() const;

  /*!
    \brief Push the initial state into the Earley parser.
    \param state The initial state to be pushed.
    \param need_expand If true, the initial state will be expanded.
  */
  void PushStateAndExpand(const ParserState& state, const bool need_expand = true);

  /*!
    \brief Reset the parser.
    \note This function is used to reset the parser, and initialize the
    parser with the root rule.
  */
  void Reset();
};

}  // namespace xgrammar

#endif  // XGRAMMAR_EARLEY_PARSER_H_
