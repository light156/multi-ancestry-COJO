#pragma once

#include <vector>
#include <cstdint>
#include <algorithm>


class LDPacked 
{
public:
    LDPacked() : n_(0) {} // default constructor

    void resize(uint64_t new_n) {
        n_ = new_n;
        data_.assign(new_n * (new_n + 1) / 2, 0.0);
    }

    inline double& operator()(uint64_t i, uint64_t j) {
        if (i > j) std::swap(i, j);
        uint64_t idx = i*n_ - i*(i-1)/2 + (j - i);
        return data_[idx];
    }

    inline double operator()(uint64_t i, uint64_t j) const {
        if (i > j) std::swap(i, j);
        uint64_t idx = i*n_ - i*(i-1)/2 + (j - i);
        return data_[idx];
    }

    inline uint64_t size() const { return n_; }

    uint64_t n_;
    std::vector<double> data_;
};
