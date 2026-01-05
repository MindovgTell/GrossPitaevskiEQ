#pragma once

// #include "ReadSimParams.hpp"
// #include "visualitsation.hpp"
// #include "SaveSimResults.hpp"


namespace gpes {

    class NonCopyable {
    protected:
        NonCopyable() = default;
        ~NonCopyable() = default;

        NonCopyable(const NonCopyable&) = delete;
        NonCopyable& operator=(const NonCopyable&) = delete;

        NonCopyable(NonCopyable&&) = delete;
        NonCopyable& operator=(NonCopyable&&) = delete;
    }; 

}