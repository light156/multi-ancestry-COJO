#pragma once
#include <vector>
#include <cstdint>
// #include <stdexcept>


class BitArray 
{
public:
    BitArray() : n_bits_(0) {}  // default constructor

    void resize(uint64_t n_bits) {
        uint64_t bits_ = (n_bits + 63) / 64;
        data_.assign(bits_, 0ULL);
        n_bits_ = n_bits;
    }

    inline void set(uint64_t idx, bool value) {
        // if (idx >= n_bits_) throw std::out_of_range("BitArray::set");
        uint64_t word = idx >> 6;
        uint64_t bit  = idx & 63;
        uint64_t mask = 1ULL << bit;
        if (value)
            data_[word] |= mask;
        else
            data_[word] &= ~mask;
    }

    inline bool get(uint64_t idx) const {
        // if (idx >= n_bits_) throw std::out_of_range("BitArray::get");
        uint64_t word = idx >> 6;
        uint64_t bit  = idx & 63;
        return (data_[word] >> bit) & 1ULL;
    }

    inline uint64_t size() const { return n_bits_; }

    std::vector<uint64_t> data_;
    uint64_t n_bits_{0};
};
