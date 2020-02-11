#include <cstdint>              // for std::int32_t
#include <cstdio>               // for std::fclose, std::fopen, std::fwrite
#include <cstdlib>              // for EXIT_FAILURE, EXIT_SUCCESS, std::exit
#include <memory>               // for std::unique_ptr
#include <numeric>              // for std::iota
#include <random>               // for std::mt19937, std::random_device
#include <string>               // for std::stoi
#include <vector>               // for std::vector
#include <boost/assert.hpp>     // for BOOST_ASSERT
#include <boost/format.hpp>     // for boost::format

namespace {
    //! A enumerated type
    /*!
    	パフォーマンスをチェックする際の対象データの種類を定義した列挙型
    */
    enum class Checktype : std::int32_t {
    	// 完全にランダムなデータ
    	RANDOM = 0,
    
    	// あらかじめソートされたデータ
    	SORT = 1,
    
    	// 最初の1/4だけソートされたデータ
    	QUARTERSORT = 2
    };
    
    //! A function.
    /*!
        ソート対象のデータを生成する
        \param checktype 生成するデータの種類
        \param n 生成するデータの個数
        \return 生成が成功したかどうか
    */
    bool make_sortdata(Checktype checktype, std::int32_t n);
}

int main(int argc, char * argv[])
{
    if (argc != 3) {
        std::exit(EXIT_FAILURE);
    }

    auto const type = static_cast<Checktype>(std::stoi(argv[1]));
    auto const n = std::stoi(argv[2]);
    
    switch (type) {
    case Checktype::RANDOM:
    	{
            if (!make_sortdata(Checktype::RANDOM, n)) {
                return EXIT_FAILURE;
            }
    	}
    	break;
    
    case Checktype::SORT:
        {
            if (!make_sortdata(Checktype::SORT, n)) {
                return EXIT_FAILURE;
            }
        }
        break;
    
    case Checktype::QUARTERSORT:
        {
            if (!make_sortdata(Checktype::QUARTERSORT, n)) {
                return EXIT_FAILURE;
            }
    	}
        break;

    default:
        BOOST_ASSERT(!"switch文のdefaultに来てしまった！");
        break;
    }

    return EXIT_SUCCESS;
}

namespace {
    bool make_sortdata(Checktype checktype, std::int32_t n)
    {
        std::vector<std::int32_t> vec(n);
        std::iota(vec.begin(), vec.end(), 1);

        std::unique_ptr< FILE, decltype(&std::fclose) > fp(nullptr, std::fclose);

        switch (checktype) {
        case Checktype::RANDOM:
        {
            std::random_device rnd;
            std::shuffle(vec.begin(), vec.end(), std::mt19937(rnd()));
            auto const filename = (boost::format("sortdata_%d_rand.dat") % n).str();
            fp = std::unique_ptr< FILE, decltype(&std::fclose) >(std::fopen(filename.c_str(), "wb"), std::fclose);
        }
        break;

        case Checktype::SORT:
        {
            auto const filename = (boost::format("sortdata_%d_already.dat") % n).str();
            fp = std::unique_ptr< FILE, decltype(&std::fclose) >(std::fopen(filename.c_str(), "wb"), std::fclose);
        }
        break;

        case Checktype::QUARTERSORT:
        {
            std::random_device rnd;
            std::shuffle(vec.begin() + n / 4, vec.end(), std::mt19937(rnd()));
            auto const filename = (boost::format("sortdata_%d_quartersort.dat") % n).str();
            fp = std::unique_ptr< FILE, decltype(&std::fclose) >(std::fopen(filename.c_str(), "wb"), std::fclose);
        }
        break;

        default:
            BOOST_ASSERT(!"switchのdefaultに来てしまった！");
            break;
        }

        if (!fp) {
            return false;
        }

        return vec.size() == std::fwrite(vec.data(), sizeof(std::int32_t), vec.size(), fp.get());
    }
}

