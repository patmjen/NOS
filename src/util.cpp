#include "util.h"

size_t std::hash<std::pair<int, int>>::operator()(const std::pair<int, int>& key) const noexcept
{
    static std::hash<long long> hasher;
    long long tohash = static_cast<long long>(key.first);
    tohash |= (static_cast<long long>(key.second) << 32);
    return hasher(tohash);
}

void graphErrFunc(const char *msg)
{
	// To avoid crashing the program we throw an exception here since Graph will otherwise call exit(1)
	throw std::runtime_error(msg);
}
