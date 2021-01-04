/*! \file parallelquicksort.cpp
    \brief スレッド並列化したクイックソートのパフォーマンスをチェックする

    Copyright © 2017-2020 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include <algorithm>                // for std::partition, std::sort
#include <array>                    // for std::array
#include <chrono>                   // for std::chrono
#include <cstdint>                  // for std::int32_t
#include <cstdio>                   // for std::fclose, std::fopen, std::fread, std::rewind
#include <execution>                // for std::execution
#include <fstream>                  // for std::ofstream
#include <iostream>                 // for std::cerr, std::cout, std::endl
#include <iterator>                 // for std::distance
#include <numeric>                  // for std::iota
#include <stack>                    // for std::stack
#include <stdexcept>                // for std::runtime_error
#include <thread>                   // for std::thread
#include <utility>                  // for std::pair
#include <vector>                   // for std::vector

#if defined(_MSC_VER) && defined(__llvm__)
    #include <string_view>          // for std::string_view_literals
    #include <system_error>         // for std::system_error 

    #define WIN32_LEAN_AND_MEAN
    #include <Windows.h>            // for GetLastError, MultiByteToWideChar
#endif

#ifndef _MSC_VER
    #include <parallel/algorithm>   // for __gnu_parallel::sort
#else
    #include <ppl.h>                // for concurrency::parallel_sort, concurrency::parallel_buffered_sort
#endif

#include <pstl/algorithm>
#include <pstl/execution>           // for pstl::execution

#include <boost/assert.hpp>         // for boost::assert
#include <boost/format.hpp>         // for boost::format
#include <boost/process.hpp>        // for boost::process
#include <boost/thread.hpp>         // for boost::thread

#include <tbb/parallel_invoke.h>    // for tbb::parallel_invoke
#include <tbb/parallel_sort.h>      // for tbb::parallel_sort

namespace {
    //! A enumerated type
    /*!
        パフォーマンスをチェックする際の対象配列の種類を定義した列挙型
    */
    enum class Checktype : std::int32_t {
        // 完全にランダムなデータ
        RANDOM = 0,

        // あらかじめソートされたデータ
        SORT = 1,

        // 最初の1/4だけソートされたデータ
        QUARTERSORT = 2
    };

    //! A global variable (constant expression).
    /*!
        計測する回数
    */
    static auto constexpr CHECKLOOP = 10;

    //! A global variable (constant expression).
    /*!
        計測する回数
    */
    static auto constexpr CPPTHREADRECMAX = 7;

    //! A global variable (constant expression).
    /*!
        ソートする配列の要素数の最初の数
    */
    static auto constexpr N = 500;

    //! A global variable (constant).
    /*!
        実行するCPUの物理コア数
    */
    static std::int32_t const NUMPHYSICALCORE = boost::thread::physical_concurrency();

    //! A global variable (constant expression).
    /*!
        クイックソートをstd::threadで並列化する際の再帰数の上限
    */
    static auto constexpr STDTHREADRECMAX = 7;

    //! A global variable (constant).
    /*!
        クイックソートにおける閾値
    */
    static auto constexpr THRESHOLD = 500;

    //! A function.
    /*!
        並列化されたソート関数のパフォーマンスをチェックする
        \param checktype パフォーマンスをチェックする際の対象配列の種類
        \param ofs 出力用のファイルストリーム
		\return 成功したかどうか
    */
    bool check_performance(Checktype checktype, std::ofstream & ofs);

    //! A function.
    /*!
        引数で与えられたstd::functionの実行時間をファイルに出力する
        \param checktype パフォーマンスをチェックする際の対象配列の種類
        \param func 実行するstd::function
        \param n 配列のサイズ
        \param ofs 出力用のファイルストリーム
        \return funcの実行結果のstd::vector
    */
    std::vector<std::int32_t> elapsed_time(Checktype checktype, std::function<void(std::vector<std::int32_t> &)> const & func, std::int32_t n, std::ofstream & ofs);

#if defined(_MSC_VER) && defined(__llvm__)
    //! A function.
    /*!
        UTF-8エンコーディングの文字列をShift-JISエンコーディングの文字列に変換する
        \param u8str UTF-8エンコーディングの文字列
        \return Shift-JISエンコーディングの文字列
    */
    std::string myutf8tosjis(std::string_view const & u8str);
#endif

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする
        \param first 範囲の下限
        \param last 範囲の上限
    */
    void quick_sort(RandomIter first, RandomIter last)
    {
        using mypair = std::pair< RandomIter, RandomIter >;

        // 範囲の情報を格納するスタック
        std::stack< mypair, std::vector< mypair > > stack;

        // 範囲の上限と下限をスタックへ積む
        stack.push(std::make_pair(first, last - 1));

        // スタックが空になるまで繰り返す
        while (!stack.empty()) {
            // 範囲の情報をスタックから取り出す
            // C++17の構造化束縛を使う
            auto const [left, right] = stack.top();
            stack.pop();

            auto i = left;
            auto j = right;
            auto const pivot = (*left + *right) / 2;

            // 左右から進めてきたiとjがぶつかるまでループ
            while (i <= j) {
                // 基準値以上の値が見つかるまで右方向へ進めていく
                while (*i < pivot) {
                    ++i;
                }

                // 基準値以下の値が見つかるまで左方向へ進めていく 
                while (*j > pivot) {
                    --j;
                }

                // 左右から進めてきたiとjがぶつかったら
                if (i <= j) {
                    // 基準値位置よりも左側にあり、基準値よりも大きい値と、
                    // 基準値位置よりも右側にあり、基準値よりも小さい値の
                    // 位置関係を交換する。
                    std::iter_swap(i, j);

                    // 次回に備えて、注目点をずらす
                    ++i;

                    // 境界チェック
                    if (j != first) {
                        // 次回に備えて、注目点をずらす
                        --j;
                    }
                }
            }

            // 左右の配列のうち、要素数が少ない方を先に処理する
            // 後で処理する側は、その範囲をスタックへ積んでおく
            if (left < j) {
                // 右側の配列を次のソート処理の対象にする
                stack.push(std::make_pair(left, j));
            }

            if (i < right) {
                // 左側の配列を次のソート処理の対象にする
                stack.push(std::make_pair(i, right));
            }
        }
    }

#if _OPENMP >= 200805
    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（OpenMPで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
        \param reci 現在の再帰の深さ
    */
    void quick_sort_openmp(RandomIter first, RandomIter last, std::int32_t reci)
    {
        // 部分ソートの要素数
        auto const num = std::distance(first, last);

        if (num <= 1) {
            // 部分ソートの要素数が1個以下なら何もすることはない
            return;
        }

        // 再帰の深さ + 1
        reci++;

        // 部分ソートが小さくなりすぎるとシリアル実行のほうが効率が良くなるため
        // 部分ソートの要素数が閾値以上の時だけ再帰させる
        // かつ、現在の再帰の深さが物理コア数以下のときだけ再帰させる
        if (num >= THRESHOLD && reci <= NUMPHYSICALCORE) {
            // 交点まで左右から入れ替えして交点を探す
            auto const middle = std::partition(first + 1, last, [first](auto n) { return n < *first; });

            // 交点 - 1の位置
            auto const mid = middle - 1;

            // 交点を移動
            std::iter_swap(first, mid);

            // 次の関数をタスクとして実行
#pragma omp task
            // 下部をソート
            quick_sort_openmp(first, mid, reci);

            // 次の関数をタスクとして実行
#pragma omp task
            // 上部をソート
            quick_sort_openmp(middle, last, reci);

            // 二つのタスクの終了を待機
#pragma omp taskwait
        }
        else {
            // 再帰なしのクイックソートの関数を呼び出す
            quick_sort(first, last);
        }
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（OpenMPで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
    */
    inline void quick_sort_openmp(RandomIter first, RandomIter last)
    {
#pragma omp parallel    // OpenMP並列領域の始まり
#pragma omp single      // task句はsingle領域で実行
        // 再帰ありの並列クイックソートを呼び出す
        quick_sort_openmp(first, last, 0);
    }
#endif

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（TBBで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
        \param reci 現在の再帰の深さ
    */
    void quick_sort_tbb(RandomIter first, RandomIter last, std::int32_t reci)
    {
        // 部分ソートの要素数
        auto const num = std::distance(first, last);

        if (num <= 1) {
            // 部分ソートの要素数が1個以下なら何もすることはない
            return;
        }

        // 再帰の深さ + 1
        reci++;

        // 部分ソートが小さくなりすぎるとシリアル実行のほうが効率が良くなるため
        // 部分ソートの要素数が閾値以上の時だけ再帰させる
        // かつ、現在の再帰の深さが物理コア数以下のときだけ再帰させる
        if (num >= THRESHOLD && reci <= NUMPHYSICALCORE) {
            // 交点まで左右から入れ替えして交点を探す
            auto const middle = std::partition(first + 1, last, [first](auto n) { return n < *first; });

            // 交点 - 1の位置
            auto const mid = middle - 1;

            // 交点を移動
            std::iter_swap(first, mid);

            // 二つのラムダ式を別スレッドで実行
            tbb::parallel_invoke(
                // 下部をソート
                [first, mid, reci]() { quick_sort_tbb(first, mid, reci); },
                // 上部をソート
                [middle, last, reci]() { quick_sort_tbb(middle, last, reci); });
        }
        else {
            // 再帰なしのクイックソートの関数を呼び出す
            quick_sort(first, last);
        }
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（TBBで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
    */
    inline void quick_sort_tbb(RandomIter first, RandomIter last)
    {
        // 再帰ありの並列クイックソートの関数を呼び出す
        quick_sort_tbb(first, last, 0);
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（std::threadで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
        \param reci 現在の再帰の深さ
    */
    void quick_sort_thread(RandomIter first, RandomIter last, std::int32_t reci)
    {
        // 部分ソートの要素数
        auto const num = std::distance(first, last);

        if (num <= 1) {
            // 部分ソートの要素数が1個以下なら何もすることはない
            return;
        }

        // 再帰の深さ + 1
        reci++;

        // 部分ソートが小さくなりすぎるとシリアル実行のほうが効率が良くなるため
        // 部分ソートの要素数が閾値以上の時だけ再帰させる
        // かつ、現在の再帰の深さが物理コア数以下のときだけ再帰させる
        if (num >= THRESHOLD && reci <= STDTHREADRECMAX) {
            // 交点まで左右から入れ替えして交点を探す
            auto const middle = std::partition(first + 1, last, [first](auto n) { return n < *first; });

            // 交点 - 1の位置
            auto const mid = middle - 1;

            // 交点を移動
            std::iter_swap(first, mid);

            // 下部をソート（別スレッドで実行）
            auto th1 = std::thread([first, mid, reci]() { quick_sort_thread(first, mid, reci); });

            // 上部をソート（別スレッドで実行）
            auto th2 = std::thread([middle, last, reci]() { quick_sort_thread(middle, last, reci); });

            // 二つのスレッドの終了を待機
            th1.join();
            th2.join();
        }
        else {
            // 再帰なしのクイックソートの関数を呼び出す
            quick_sort(first, last);
        }
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（std::threadで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
    */
    inline void quick_sort_thread(RandomIter first, RandomIter last)
    {
        // 再帰ありの並列クイックソートの関数を呼び出す
        quick_sort_thread(first, last, 0);
    }

#ifdef DEBUG_CHECK_RECUL
    //! A function.
    /*!
        再帰数の上限でパフォーマンスがどう変化するかチェックする
        \param ofs 出力用のファイルストリーム
    */
    void check_performance_recul(std::ofstream & ofs);

    //! A function.
    /*!
        引数で与えられたstd::functionの実行時間をファイルに出力する
        \param func 実行するstd::function
        \param n 配列のサイズ
        \param ofs 出力用のファイルストリーム
    */
    void elapsed_time_recul(std::function<void(std::vector<std::int32_t> &)> const & func, std::int32_t n, std::ofstream & ofs);

#if _OPENMP >= 200805
    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（OpenMPで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
        \param reci 現在の再帰の深さ
        \param recul 再帰数の上限
    */
    void quick_sort_openmp_recul(RandomIter first, RandomIter last, std::int32_t reci, std::int32_t recul)
    {
        // 部分ソートの要素数
        auto const num = std::distance(first, last);

        if (num <= 1) {
            // 部分ソートの要素数が1個以下なら何もすることはない
            return;
        }

        // 再帰の深さ + 1
        reci++;

        // 部分ソートが小さくなりすぎるとシリアル実行のほうが効率が良くなるため
        // 部分ソートの要素数が閾値以上の時だけ再帰させる
        // かつ、現在の再帰の深さが物理コア数以下のときだけ再帰させる
        if (num >= THRESHOLD && reci <= recul) {
            // 交点まで左右から入れ替えして交点を探す
            auto const middle = std::partition(first + 1, last, [first](auto n) { return n < *first; });

            // 交点 - 1の位置
            auto const mid = middle - 1;

            // 交点を移動
            std::iter_swap(first, mid);

            // 次の関数をタスクとして実行
#pragma omp task
            // 下部をソート
            quick_sort_openmp_recul(first, mid, reci, recul);

            // 次の関数をタスクとして実行
#pragma omp task
            // 上部をソート
            quick_sort_openmp_recul(middle, last, reci, recul);

            // 二つのタスクの終了を待機
#pragma omp taskwait
        }
        else {
            // 再帰なしのクイックソートの関数を呼び出す
            quick_sort(first, last);
        }
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（OpenMPで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
        \param recul 再帰数の上限
    */
    inline void quick_sort_openmp_recul(RandomIter first, RandomIter last, std::int32_t recul)
    {
#pragma omp parallel    // OpenMP並列領域の始まり
#pragma omp single      // task句はsingle領域で実行
        // 再帰ありの並列クイックソートを呼び出す
        quick_sort_openmp_recul(first, last, 0, recul);
    }
#endif

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（TBBで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
        \param reci 現在の再帰の深さ
        \param recul 再帰数の上限
    */
    void quick_sort_tbb_recul(RandomIter first, RandomIter last, std::int32_t reci, std::int32_t recul)
    {
        // 部分ソートの要素数
        auto const num = std::distance(first, last);

        if (num <= 1) {
            // 部分ソートの要素数が1個以下なら何もすることはない
            return;
        }

        // 再帰の深さ + 1
        reci++;

        // 部分ソートが小さくなりすぎるとシリアル実行のほうが効率が良くなるため
        // 部分ソートの要素数が閾値以上の時だけ再帰させる
        // かつ、現在の再帰の深さが物理コア数以下のときだけ再帰させる
        if (num >= THRESHOLD && reci <= recul) {
            // 交点まで左右から入れ替えして交点を探す
            auto const middle = std::partition(first + 1, last, [first](auto n) { return n < *first; });

            // 交点 - 1の位置
            auto const mid = middle - 1;

            // 交点を移動
            std::iter_swap(first, mid);

            // 二つのラムダ式を別スレッドで実行
            tbb::parallel_invoke(
                // 下部をソート
                [first, mid, reci, recul]() { quick_sort_tbb_recul(first, mid, reci, recul); },
                // 上部をソート
                [middle, last, reci, recul]() { quick_sort_tbb_recul(middle, last, reci, recul); });
        }
        else {
            // 再帰なしのクイックソートの関数を呼び出す
            quick_sort(first, last);
        }
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（TBBで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
        \param recul 再帰数の上限
    */
    inline void quick_sort_tbb_recul(RandomIter first, RandomIter last, std::int32_t recul)
    {
        // 再帰ありの並列クイックソートの関数を呼び出す
        quick_sort_tbb_recul(first, last, 0, recul);
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（std::threadで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
        \param reci 現在の再帰の深さ
        \param recul 再帰数の上限
    */
    void quick_sort_thread_recul(RandomIter first, RandomIter last, std::int32_t reci, std::int32_t recul)
    {
        // 部分ソートの要素数
        auto const num = std::distance(first, last);

        if (num <= 1) {
            // 部分ソートの要素数が1個以下なら何もすることはない
            return;
        }

        // 再帰の深さ + 1
        reci++;

        // 部分ソートが小さくなりすぎるとシリアル実行のほうが効率が良くなるため
        // 部分ソートの要素数が閾値以上の時だけ再帰させる
        // かつ、現在の再帰の深さが物理コア数以下のときだけ再帰させる
        if (num >= THRESHOLD && reci <= recul) {
            // 交点まで左右から入れ替えして交点を探す
            auto const middle = std::partition(first + 1, last, [first](auto n) { return n < *first; });

            // 交点 - 1の位置
            auto const mid = middle - 1;

            // 交点を移動
            std::iter_swap(first, mid);

            // 下部をソート（別スレッドで実行）
            auto th1 = std::thread([first, mid, reci, recul]() { quick_sort_thread_recul(first, mid, reci, recul); });

            // 上部をソート（別スレッドで実行）
            auto th2 = std::thread([middle, last, reci, recul]() { quick_sort_thread_recul(middle, last, reci, recul); });

            // 二つのスレッドの終了を待機
            th1.join();
            th2.join();
        }
        else {
            // 再帰なしのクイックソートの関数を呼び出す
            quick_sort(first, last);
        }
    }

    template < class RandomIter >
    //! A template function.
    /*!
        指定された範囲の要素をクイックソートでソートする（std::threadで並列化）
        \param first 範囲の下限
        \param last 範囲の上限
    */
    inline void quick_sort_thread_recul(RandomIter first, RandomIter last, std::int32_t recul)
    {
        // 再帰ありの並列クイックソートの関数を呼び出す
        quick_sort_thread_recul(first, last, 0, recul);
    }
#endif

#ifdef DEBUG
    //! A template function.
    /*!
        与えられた二つのstd::vectorのすべての要素が同じかどうかチェックする
        \param v1 一つ目のstd::vector
        \param v2 二つめのstd::vector
        \return 与えられた二つのstd::vectorのすべての要素が同じならtrue、そうでないならfalse
    */
    bool vec_check(std::vector<std::int32_t> const & v1, std::vector<std::int32_t> const & v2);
#endif
}

int main()
{
#if defined(_MSC_VER) && defined(__llvm__)
    using namespace std::string_view_literals;

    std::cout << myutf8tosjis(u8R"(物理コア数: )"sv) << boost::thread::physical_concurrency();
    std::cout << myutf8tosjis(u8R"(, 論理コア数: )"sv) << boost::thread::hardware_concurrency() << std::endl;

#ifdef DEBUG_CHECK_RECUL
    std::ofstream ofsrec(myutf8tosjis(u8R"(完全にシャッフルされたデータ_再帰数チェック.csv)"sv));
    check_performance_recul(ofsrec);
#else
    std::ofstream ofsrandom(myutf8tosjis(u8R"(完全にシャッフルされたデータ.csv)"sv));
    std::ofstream ofssort(myutf8tosjis(u8R"(あらかじめソートされたデータ.csv)"sv));
    std::ofstream ofsquartersort(myutf8tosjis(u8R"(最初の1_4だけソートされたデータ.csv)"sv));

    std::cout << myutf8tosjis(u8R"(完全にシャッフルされたデータを計測中...)"sv) << '\n';
    if (!check_performance(Checktype::RANDOM, ofsrandom)) {
        return -1;
    }

    std::cout << '\n' << myutf8tosjis(u8R"(あらかじめソートされたデータを計測中...)"sv) << '\n';
    if (!check_performance(Checktype::SORT, ofssort)) {
        return -1;
    }

    std::cout << '\n' << myutf8tosjis(u8R"(最初の1/4だけソートされたデータを計測中...)"sv) << '\n';
    if (!check_performance(Checktype::QUARTERSORT, ofsquartersort)) {
        return -1;
    }
#endif

#else
    std::cout << "物理コア数: " << boost::thread::physical_concurrency();
    std::cout << ", 論理コア数: " << boost::thread::hardware_concurrency() << std::endl;

#ifdef DEBUG_CHECK_RECUL
    std::ofstream ofsrec("完全にシャッフルされたデータ_再帰数チェック.csv");
    check_performance_recul(ofsrec);
#else
    std::ofstream ofsrandom("完全にシャッフルされたデータ.csv");
    std::ofstream ofssort("あらかじめソートされたデータ.csv");
    std::ofstream ofsquartersort("最初の1_4だけソートされたデータ.csv");

    std::cout << "完全にシャッフルされたデータを計測中...\n";
    if (!check_performance(Checktype::RANDOM, ofsrandom)) {
        return -1;
    }

    std::cout << "\nあらかじめソートされたデータを計測中...\n";
    if (!check_performance(Checktype::SORT, ofssort)) {
        return -1;
    }

    std::cout << "\n最初の1/4だけソートされたデータを計測中...\n";
    if (!check_performance(Checktype::QUARTERSORT, ofsquartersort)) {
        return -1;
    }
#endif
#endif

    return 0;
}

namespace {
    bool check_performance(Checktype checktype, std::ofstream & ofs)
    {
#ifndef _MSC_VER
        std::array< std::uint8_t, 3 > const bom = { 0xEF, 0xBB, 0xBF };
        ofs.write(reinterpret_cast<char const *>(bom.data()), sizeof(bom));
#endif

#if defined(_MSC_VER) && _OPENMP < 200805
        ofs << "配列の要素数,std::sort,クイックソート,std::thread,TBB,concurrency::parallel_sort,concurrency::parallel_buffered_sort,tbb::parallel_sort,std::sort (MSVC内蔵のParallelism TS),std::sort (Parallel STLのParallelism TS)\n";
#elif defined(_MSC_VER)
        ofs << "配列の要素数,std::sort,クイックソート,std::thread,TBB,OpenMP,concurrency::parallel_sort,concurrency::parallel_buffered_sort,tbb::parallel_sort,std::sort (MSVC内蔵のParallelism TS),std::sort (Parallel STLのParallelism TS)\n";
#else
        ofs << u8R"(配列の要素数,std::sort,クイックソート,std::thread,OpenMP,TBB,__gnu_parallel::sort,tbb::parallel_sort,std::sort (Parallelism TS),std::sort (Parallel STLのParallelism TS))" << '\n';
#endif
        
        auto issuccess = true;

        auto n = N;
        for (auto i = 0; i < 6; i++) {
            for (auto j = 0; j < 2; j++) {
#if defined(_MSC_VER) && defined(__llvm__)
                using namespace std::string_view_literals;

                std::cout << n << myutf8tosjis(u8R"(個を計測中...)"sv) << '\n';
#else
                std::cout << n << "個を計測中...\n";
#endif
                std::vector<std::int32_t> vec(n);
                std::iota(vec.begin(), vec.end(), 1);

                std::array< std::vector<std::int32_t>, 10 > vecar;

                ofs << n << ',';

                try {
                    vecar[0] = elapsed_time(checktype, [](auto && vec) { std::sort(vec.begin(), vec.end()); }, n, ofs);
                    vecar[1] = elapsed_time(checktype, [](auto && vec) { quick_sort(vec.begin(), vec.end()); }, n, ofs);
                    vecar[2] = elapsed_time(checktype, [](auto && vec) { quick_sort_thread(vec.begin(), vec.end()); }, n, ofs);
                
#if _OPENMP >= 200805
                    vecar[3] = elapsed_time(checktype, [](auto && vec) { quick_sort_openmp(vec.begin(), vec.end()); }, n, ofs);
#endif
                    vecar[4] = elapsed_time(checktype, [](auto && vec) { quick_sort_tbb(vec.begin(), vec.end()); }, n, ofs);

#ifndef _MSC_VER
                    vecar[5] = elapsed_time(checktype, [](auto && vec) { __gnu_parallel::sort(vec.begin(), vec.end()); }, n, ofs);
#else
                    vecar[5] = elapsed_time(checktype, [](auto && vec) { concurrency::parallel_sort(vec.begin(), vec.end()); }, n, ofs);

                    vecar[6] = elapsed_time(checktype, [](auto && vec) { concurrency::parallel_buffered_sort(vec.begin(), vec.end()); }, n, ofs);
#endif


                    vecar[7] = elapsed_time(checktype, [](auto && vec) { tbb::parallel_sort(vec); }, n, ofs);

                    vecar[8] = elapsed_time(checktype, [](auto && vec) { std::sort(std::execution::par, vec.begin(), vec.end()); }, n, ofs);

#ifdef _MSC_VER
                    vecar[9] = elapsed_time(checktype, [](auto && vec) { std::sort(pstl::execution::par, vec.begin(), vec.end()); }, n, ofs);
#else
                    vecar[9] = elapsed_time(checktype, [](auto && vec) { std::sort(__pstl::execution::par, vec.begin(), vec.end()); }, n, ofs);
#endif

                }
                catch (std::runtime_error const & e) {
                    std::cerr << e.what() << std::endl;
                    return false;
                }

                ofs << std::endl;

#ifdef DEBUG
                for (auto k = 0U; k < vecar.size(); k++) {
#if _OPENMP < 200805
                    if (k == 3) {
                        continue;
                    }
#endif

                    if (!k || k == 6 || k == 7 || k == 8) {
                        continue;
                    }

                    if (static_cast<std::int32_t>(vecar[k].size()) != n) {
                        issuccess = false;
                        continue;
                    }

                    if (!vec_check(vec, vecar[k])) {
                        std::cerr << "Error! vecar[" << k << ']' << std::endl;
                        issuccess = false;
                    }
                }
#endif
                if (!j) {
                    n *= 2;
                }
            }

            n *= 5;
        }

        return issuccess;
    }

    std::vector<std::int32_t> elapsed_time(Checktype checktype, std::function<void(std::vector<std::int32_t> &)> const & func, std::int32_t n, std::ofstream & ofs)
    {
        using namespace std::chrono;

        std::vector<std::int32_t> vec(n);
        std::unique_ptr< FILE, decltype(&std::fclose) > fp(nullptr, fclose);

        auto const program_name = "makequicksortdata";

        switch (checktype) {
        case Checktype::RANDOM:
            {
                auto const filename = (boost::format("sortdata_%d_rand.dat") % n).str();
                fp = std::unique_ptr< FILE, decltype(&std::fclose) >(std::fopen(filename.c_str(), "rb"), std::fclose);
                if (!fp) {
                    boost::process::child(program_name + (boost::format(" 0 %d") % n).str()).wait();
                    fp = std::unique_ptr< FILE, decltype(&std::fclose) >(std::fopen(filename.c_str(), "rb"), std::fclose);
                }
            }
            break;

        case Checktype::SORT:
            {
                auto const filename = (boost::format("sortdata_%d_already.dat") % n).str();
                fp = std::unique_ptr< FILE, decltype(&std::fclose) >(std::fopen(filename.c_str(), "rb"), std::fclose);
                if (!fp) {
                    boost::process::child(program_name + (boost::format(" 1 %d") % n).str()).wait();
                    fp = std::unique_ptr< FILE, decltype(&std::fclose) >(std::fopen(filename.c_str(), "rb"), std::fclose);
                }
            }
            break;

        case Checktype::QUARTERSORT:
            {
                auto const filename = (boost::format("sortdata_%d_quartersort.dat") % n).str();
                fp = std::unique_ptr< FILE, decltype(&std::fclose) >(std::fopen(filename.c_str(), "rb"), std::fclose);
                if (!fp) {
                    boost::process::child(program_name + (boost::format(" 2 %d") % n).str()).wait();
                    fp = std::unique_ptr< FILE, decltype(&std::fclose) >(std::fopen(filename.c_str(), "rb"), std::fclose);
                }
            }
            break;

        default:
            BOOST_ASSERT(!"switch文のdefaultに来てしまった！");
            break;
        }

        auto const readsize = vec.size();        
        if (readsize != std::fread(vec.data(), sizeof(std::int32_t), readsize, fp.get())) {
            throw std::runtime_error("std::freadに失敗");
        }

        auto elapsed_time = 0.0;
        
        for (auto i = 1; i <= CHECKLOOP; i++) {
            auto const beg = high_resolution_clock::now();
            func(vec);
            auto const end = high_resolution_clock::now();

            elapsed_time += (duration_cast<duration<double>>(end - beg)).count();
            
            if (i != CHECKLOOP) {
                std::rewind(fp.get());
                if (readsize != std::fread(vec.data(), sizeof(std::int32_t), readsize, fp.get())) {
                    throw std::runtime_error("std::freadに失敗");
                }
            }
        }

        ofs << boost::format("%.10f") % (elapsed_time / static_cast<double>(CHECKLOOP)) << ',';

        return vec;
    }

#ifdef DEBUG_CHECK_RECUL
    void check_performance_recul(std::ofstream& ofs)
    {
        using namespace std::string_view_literals;

#if _OPENMP < 200805
        ofs << "再帰数,std::thread,TBB\n";
#else
        ofs << myutf8tosjis(u8R"(再帰数,std::thread,OpenMP,TBB)"sv) << '\n';
#endif

        for (auto recul = 0U; recul <= boost::thread::hardware_concurrency(); recul++) {
#if defined(_MSC_VER) && defined(__llvm__)
            std::cout << myutf8tosjis(u8R"(再帰数: )"sv) << recul << myutf8tosjis(u8R"(を計測中...)"sv) << '\n';
#else
            std::cout << "再帰数: " << recul << "を計測中...\n";
#endif

            ofs << recul << ',';
            elapsed_time_recul([recul](auto && vec) { quick_sort_thread_recul(vec.begin(), vec.end(), recul); }, 5000000, ofs);
#if _OPENMP >= 200805
            elapsed_time_recul([recul](auto && vec) { quick_sort_openmp_recul(vec.begin(), vec.end(), recul); }, 5000000, ofs);
#endif
            elapsed_time_recul([recul](auto && vec) { quick_sort_tbb_recul(vec.begin(), vec.end(), recul); }, 5000000, ofs);
            ofs << std::endl;
        }
    }

    void elapsed_time_recul(std::function<void(std::vector<std::int32_t>&)> const& func, std::int32_t n, std::ofstream& ofs)
    {
        using namespace std::chrono;

        std::vector<std::int32_t> vec(n);
        std::unique_ptr< FILE, decltype(&std::fclose) > fp(nullptr, fclose);

        auto const program_name = "makequicksortdata";

        auto const filename = (boost::format("sortdata_%d_rand.dat") % n).str();
        fp = std::unique_ptr< FILE, decltype(&std::fclose) >(std::fopen(filename.c_str(), "rb"), std::fclose);
        if (!fp) {
            boost::process::child(program_name + (boost::format(" 0 %d") % n).str()).wait();
            fp = std::unique_ptr< FILE, decltype(&std::fclose) >(std::fopen(filename.c_str(), "rb"), std::fclose);
        }

        auto const readsize = vec.size();
        if (readsize != std::fread(vec.data(), sizeof(std::int32_t), readsize, fp.get())) {
            throw std::runtime_error("std::freadに失敗");
        }

        auto elapsed_time = 0.0;

        for (auto i = 1; i <= CHECKLOOP; i++) {
            auto const beg = high_resolution_clock::now();
            func(vec);
            auto const end = high_resolution_clock::now();

            elapsed_time += (duration_cast<duration<double>>(end - beg)).count();

            if (i != CHECKLOOP) {
                std::rewind(fp.get());
                if (readsize != std::fread(vec.data(), sizeof(std::int32_t), readsize, fp.get())) {
                    throw std::runtime_error("std::freadに失敗");
                }
            }
        }

        ofs << boost::format("%.10f") % (elapsed_time / static_cast<double>(CHECKLOOP)) << ',';
    }
#endif

#if defined(_MSC_VER) && defined(__llvm__)
    std::string myutf8tosjis(std::string_view const & u8str)
    {
        // UTF-8 -> UTF-16
        auto length = ::MultiByteToWideChar(
            CP_UTF8,
            0,
            reinterpret_cast<const char*>(u8str.data()),
            static_cast<int>(u8str.length()),
            nullptr,
            0);
        if (!length) {
            throw std::system_error(std::error_code(GetLastError(), std::system_category()));
        }

        std::wstring temp(length, '\0');

        auto res = ::MultiByteToWideChar(
            CP_UTF8,
            0,
            reinterpret_cast<const char*>(u8str.data()),
            static_cast<int>(u8str.length()),
            temp.data(), temp.length());
        if (!res) {
            throw std::system_error(std::error_code(GetLastError(), std::system_category()));
        }

        // UTF-16 -> Shift-JIS
        length = ::WideCharToMultiByte(CP_ACP, 0,
            temp.data(), static_cast<int>(temp.length()),
            nullptr, 0,
            nullptr, nullptr);

        std::string result(length, '\0');

        res = ::WideCharToMultiByte(CP_ACP, 0,
            temp.data(), static_cast<int>(temp.length()),
            result.data(), static_cast<int>(result.length()),
            nullptr, nullptr);
        if (!res) {
            throw std::system_error(std::error_code(GetLastError(), std::system_category()));
        }

        return result;
    }
#endif

#ifdef DEBUG
    bool vec_check(std::vector<std::int32_t> const & v1, std::vector<std::int32_t> const & v2)
    {
        auto const size = v1.size();
        BOOST_ASSERT(size == v2.size());

        for (auto i = 0UL; i < size; i++) {
            if (v1[i] != v2[i]) {
                std::cerr << "Error! i = " << i << '\n';
                return false;
            }
        }

        return true;
    }
#endif
}
