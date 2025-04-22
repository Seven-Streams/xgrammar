/*!
 *  Copyright (c) 2025 by Contributors
 * \file xgrammar/earley_parser.h
 * \brief The header for the definition of the Earley parser.
 */
#ifndef XGRAMMAR_EARLEY_PARSER_H_
#define XGRAMMAR_EARLEY_PARSER_H_
#include <cstdint>
#include <vector>

#include "xgrammar/grammar.h"
namespace xgrammar {
// Be the same as StackElement.
struct State {
  /*! \brief The rule's id. Used for debug purposes. */
  int32_t rule_id = -1;
  /*! \brief Which choice in this rule is selected. */
  int32_t sequence_id = -1;
  /*! \brief Which element of the choice sequence is to be visited. When the current sequence is
   * a tag dispatch rule, this element id the currently visited node. */
  int32_t element_id = -1;

  /*! \brief The number of left utf8 bytes in the current element. Used when the element is
   * a character class or a character class star. */
  int32_t left_utf8_bytes = 0;
  /*! \brief The next position to match in the current byte string. Used when the element is
   * a byte string. */
  int32_t element_in_string = 0;

  /*! \brief The id of the parent node in the PersistentStack. */
  int32_t parent_id = -1;

  /*! \brief The reference count of this StackElement. If reduces to zero, the node will be
   * removed from the StackElementBuffer. */
  int reference_count = 0;

  /*! \brief A parent_id value of kNoParent means this StackElement is the root of the tree. */
  static constexpr int32_t kNoParent = -1;

  constexpr State() = default;
  constexpr State(
      int32_t rule_id, int32_t sequence_id, int32_t element_id, int32_t parent_id = kNoParent
  )
      : rule_id(rule_id), sequence_id(sequence_id), element_id(element_id), parent_id(parent_id) {}

  // The element is invalid when sequence_id is -1.
  bool IsInvalid() const { return sequence_id == -1; }

  bool operator==(const State& other) const {
    return rule_id == other.rule_id && sequence_id == other.sequence_id &&
           element_id == other.element_id && parent_id == other.parent_id &&
           left_utf8_bytes == other.left_utf8_bytes && element_in_string == other.element_in_string;
  }
};
class EarleyParser {
 private:
  /*! \brief The grammar to be parsed. */
  Grammar grammar_;
  /*! \brief The tree storing all states. It's used for completation. */
  std::vector<std::vector<State>> states;
  /*!
      \brief The history of states. i.e. the i-th(0-base) vector
      will store the states after matching i characters. It's used
      for rollback.
   */
  std::vector<std::vector<State>> history_states;

  /*!
    \brief The scanning operation of the Earley parser.
  */
  bool Scan(const State& state, const uint8_t& ch) const;

  /*!
      \brief The prediction operation of the Earley parser.
  */
  void Predict(const State& state) const;

  /*!
      \brief The completion operation of the Earley parser.
  */
  void Complete(const State& state) const;

 public:
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
  void PopBackStates(int32_t count);

  /*!
    \brief Check if the root rule is completed.
    \return True if the root rule is completed, false otherwise.
  */
  bool CanReachEnd() const;
};
}  // namespace xgrammar
#endif  // XGRAMMAR_EARLEY_PARSER_H_
