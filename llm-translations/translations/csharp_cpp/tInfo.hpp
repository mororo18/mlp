#pragma once

namespace MLP {

    class tInfoIndex {
    public:
        static const int T = 0;
        static const int C = 1;
        static const int W = 2;
    };

    class Neighborhood {
    public:
        static const int SWAP        = 0;
        static const int REINSERTION = 1;
        static const int OR_OPT2     = 2;
        static const int OR_OPT3     = 3;
        static const int TWO_OPT     = 4;
    };

}
