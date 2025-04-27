#ifndef XGRAMMAR_ELEMENT_BUFFER_H_
#define XGRAMMAR_ELEMENT_BUFFER_H_
#include <cstdint>
#include <vector>

#include "logging.h"
namespace xgrammar {
template <class T, const T& invalid_element>
class ElementBuffer {
 public:
  /*!
   * \brief Allocate a new StackElement. with given initial value.
   * \returns The id of the allocated node.
   */
  int32_t Allocate(T stack_element) {
    int32_t id;
    if (free_nodes_.empty()) {
      buffer_.emplace_back();
      id = static_cast<int32_t>(buffer_.size()) - 1;
    } else {
      id = free_nodes_.back();
      XGRAMMAR_DCHECK(buffer_[id].IsInvalid());
      free_nodes_.pop_back();
    }
    stack_element.reference_count = 0;
    buffer_[id] = stack_element;
    return id;
  }

  /*! \brief Free the StackElement with the given id. */
  void Free(int32_t id) {
    XGRAMMAR_DCHECK(!buffer_[id].IsInvalid());
    buffer_[id] = invalid_element;
    free_nodes_.push_back(id);
  }

  /*! \brief Get the capacity of the buffer. */
  size_t Capacity() const { return buffer_.size(); }

  /*! \brief Get the number of allocated nodes. */
  size_t Size() const {
    XGRAMMAR_DCHECK(buffer_.size() >= free_nodes_.size());
    return buffer_.size() - free_nodes_.size();
  }

  /*! \brief Get the StackElement with the given id. */
  T& operator[](int32_t id) {
    XGRAMMAR_DCHECK(id >= 0 && id < static_cast<int32_t>(buffer_.size()));
    XGRAMMAR_DCHECK(!buffer_[id].IsInvalid());
    return buffer_[id];
  }
  const T& operator[](int32_t id) const {
    XGRAMMAR_DCHECK(id >= 0 && id < static_cast<int32_t>(buffer_.size()));
    XGRAMMAR_DCHECK(!buffer_[id].IsInvalid());
    return buffer_[id];
  }

  void Reset() {
    buffer_.clear();
    free_nodes_.clear();
  }

 private:
  /*! \brief The buffer to store all States. */
  std::vector<T> buffer_;
  /*! \brief A stack to store all free node ids. */
  std::vector<int32_t> free_nodes_;
};
}  // namespace xgrammar
#endif  // XGRAMMAR_ELEMENT_BUFFER_H_
