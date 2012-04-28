#include <iostream>
#include "common.h"
#include "bsconfig.h"

using namespace std;

long BSConfig::totalNumOfBlocks() {
    return min((M/WINDOW_WIDTH), BLOCK_SIZE);
}

long BSConfig::totalNumOfThread() {
    return totalNumOfBlocks() * WINDOW_WIDTH;
}

BSConfig::BSConfig() {
    S = 0.0; E = 0.0; r = 0.0; sigma = 0.0;
    T = 0.0;
    M = 0;
    DEBUG_LEVEL = false;
    RND_MODE = 0;
}
