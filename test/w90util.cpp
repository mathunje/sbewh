#include "gtest/gtest.h"

#include<vector>
#include "w90util.h"

namespace{

TEST(w90util, readTB){
    const char *fname = "../data/example_tb.dat";
    W90_tb wtb = w90_read_tb(fname, true);
    ASSERT_TRUE(wtb.valid);
    int foundCnt = 0;
    for(unsigned i=0; i<wtb.wsCells.size(); i++){
        if ( wtb.wsCells[i].r[0] == 0 && wtb.wsCells[i].r[1] == 0 && wtb.wsCells[i].r[2] == 0){
            foundCnt++;
            EXPECT_NEAR( std::real(wtb.wsCells[i].H(0, 1)), atomicUnits::from_eV(0.14439922E+00), 1e-10);
            EXPECT_NEAR( std::imag(wtb.wsCells[i].H(0, 1)), atomicUnits::from_eV(0.22271074E-12), 1e-10);
            EXPECT_NEAR( std::real(wtb.wsCells[i].D[1](1, 0)), atomicUnits::from_A(-0.87281205E-01), 1e-10);
            EXPECT_NEAR( std::imag(wtb.wsCells[i].D[1](1, 0)), atomicUnits::from_A(0.80161572E-14), 1e-10);
        }
    }
    EXPECT_EQ(foundCnt, 1);
}

TEST(w90util, readWsvec){
    const char *fname = "../data/example_wsvec.dat";
    W90_wsvec wsvec = w90_read_wsvec(fname, true);
    ASSERT_TRUE(wsvec.valid);
    EXPECT_EQ(wsvec.cells.size(), 21 * 2 * 2);
    unsigned ti = 5 * 4 + 1, tj = 20; // line 84
    W90_wsvecCell c = wsvec.cells[ti];
    EXPECT_EQ(c.r[0], -1);
    EXPECT_EQ(c.r[1], 0);
    EXPECT_EQ(c.r[2], 1);
    EXPECT_EQ(c.m, 0);
    EXPECT_EQ(c.n, 1);
    c = wsvec.cells[20];
    EXPECT_EQ(c.T.size(), 4);
    EXPECT_EQ(c.T[2][0], 2);
    EXPECT_EQ(c.T[2][1], 0);
    EXPECT_EQ(c.T[2][2], -2);
}

TEST(w90util, symmetrize){
    W90_tb w;
    w.wsCells.resize(3, 4);
    auto &f1 = w.wsCells[0];
    auto &f2 = w.wsCells[1];
    auto &z = w.wsCells[2];
    f1.r = {0, 0, 1};
    f2.r = -f1.r;
    z.r = {0, 0, 0};
    f1.H(0, 1) = 2;
    f2.H(1, 0) = std::complex<double>(1, -10);
    z.H(0, 0) = std::complex<double>(1, 224);
    z.H(1, 2) = std::complex<double>(1, 5);
    z.H(2, 1) = std::complex<double>(2, -3);

    for(unsigned dir=0; dir<3; dir++)
        z.D[dir](dir, dir) = dir+1;

    ASSERT_TRUE( w90_symmetrize(w, true) );
    EXPECT_EQ ( std::real(f1.H(0, 1)), 1.5);
    EXPECT_EQ ( std::imag(f2.H(1, 0)), -5);

    EXPECT_EQ ( std::real(z.H(0, 0)), 1);
    EXPECT_EQ ( std::imag(z.H(0, 0)), 0);

    EXPECT_EQ ( std::real(z.H(1, 2)), 1.5);
    EXPECT_EQ ( std::imag(z.H(2, 1)), -4);

    EXPECT_EQ( std::real(z.D[0](0, 0)), 1);
    EXPECT_EQ( std::real(z.D[1](1, 1)), 2);
    EXPECT_EQ( std::real(z.D[2](2, 2)), 3);
}

/*******************************************/

W90_tb createTestTb(unsigned numWann){
    W90_tb tb;
    tb.valid = true;
    tb.numWann = numWann;
    W90_wignerSeitzCell wsc0(numWann);
    wsc0.r = CellIndex({0, 0, 0});
    wsc0.degeneracy = 1;
    for(int i=0; i<numWann; i++)
        for(int j=0; j<numWann; j++){
            wsc0.H(i, j) = i*numWann + j;
            for(int d=0; d<3; d++)
                wsc0.D[d](i, j) = d + i*numWann +j;
        }
    tb.wsCells.push_back(wsc0);
    W90_wignerSeitzCell wsc1(numWann);
    wsc1.r = CellIndex({1, 1, 1});
    wsc1.degeneracy = 2;
    for(int i=0; i<numWann; i++)
        for(int j=0; j<numWann; j++){
            wsc1.H(i, j) = 10 + i*numWann + j;
            for(int d=0; d<3; d++)
                wsc1.D[d](i, j) = 10 + d + i*numWann +j;
        }
    tb.wsCells.push_back(wsc1);
    return tb;
}

TEST(w90util, oldRealSpaceOperators){
    unsigned numWann = 2;
    W90_tb tb = createTestTb(numWann);
    auto rso = w90_calcRealSpaceOperators(tb);
    auto it = rso.find({0, 0, 0});
    ASSERT_TRUE( it!=rso.end() );
    EXPECT_NEAR( std::real(it->second.H(1, 0)), 1*numWann + 0, 1e-10);
    EXPECT_NEAR( std::real(it->second.D[2](1, 0)), 2 + 1*numWann + 0, 1e-10);
    it = rso.find({1, 1, 1});
    ASSERT_TRUE( it!=rso.end() );
    EXPECT_NEAR( std::real(it->second.H(0, 1)), (double)(10 + 0*numWann + 1) / 2, 1e-10);
    EXPECT_NEAR( std::real(it->second.D[2](0, 1)), (double)(10 + 2 + 0*numWann + 1) / 2, 1e-10);
}

TEST(w90util, newRealSpaceOperators){
    unsigned numWann = 2;
    W90_tb tb = createTestTb(numWann);
    W90_wsvec wsvec;
    wsvec.valid = true;
    W90_wsvecCell c;
    c.r = CellIndex({0, 0, 0});
    for(unsigned m=0; m<numWann; m++)
        for(unsigned n=0; n<numWann; n++){
            c.m = m;
            c.n = n;
            c.T = std::vector<CellIndex>({ {2, 2, 2}, {0, 0, 0} });
            wsvec.cells.push_back(c);
        }
    c.r = CellIndex({1, 1, 1});
    for(unsigned m=0; m<numWann; m++)
        for(unsigned n=0; n<numWann; n++){
            c.m = m;
            c.n = n;
            c.T.clear();
            c.T.push_back({0, 0, 0});
            if ( m )
                c.T.push_back({1, 1, 1});
            wsvec.cells.push_back(c);
        }
    std::unordered_map<CellIndex, W90_realSpaceOperators> rso = w90_calcRealSpaceOperators(tb, wsvec);
    ASSERT_EQ( rso.size(), 3);
    auto it = rso.find({0, 0, 0});
    ASSERT_TRUE( it != rso.end() );
    EXPECT_NEAR( std::real(it->second.H(1, 0)), 0.5 * (1*numWann + 0), 1e-10);
    EXPECT_NEAR( std::real(it->second.D[2](1, 0)), 0.5 * (2 + 1*numWann + 0), 1e-10);
    it = rso.find({2, 2, 2});
    ASSERT_TRUE( it != rso.end() );
    EXPECT_NEAR( std::real(it->second.H(0, 1)), 0.5 * (0*numWann +1), 1e-10);
    EXPECT_NEAR( std::real(it->second.D[2](0, 1)), 0.5 * (2 + 0*numWann+1), 1e-10);

    EXPECT_NEAR( std::real(it->second.H(1, 0)), 0.25 * 10 + 0.75 *(1*numWann + 0), 1e-10);
    EXPECT_NEAR( std::real(it->second.D[2](1, 0)), 0.25 * 10 + 0.75 *(2 + 1*numWann + 0), 1e-10);
}

TEST(w90util, readUmat){
    auto [Uks, success] = w90_read_Umat("../data/CdSe_u_dis.mat", true);
    ASSERT_TRUE ( success );
    ASSERT_EQ ( Uks.size(), 27);
    EXPECT_EQ ( Uks[0].U.getRowSize(), 12);
    EXPECT_EQ ( Uks[0].U.getColSize(), 8);
    for(unsigned d=0; d<3; d++)
        EXPECT_NEAR(Uks[0].k[d], -1.0/3, 1e-6);
    for(unsigned s=0; s<Uks[0].U.getColSize(); s++)
        EXPECT_FLOAT_EQ( std::real(Uks[0].U(s, s)), 1);
    std::complex<double> entry{-0.0401003833, -0.0058553036};
    std::complex<double> readValue = Uks[25].U(Uks[25].U.getRowSize()-1, Uks[25].U.getColSize()-1);
    EXPECT_FLOAT_EQ( std::real(entry), std::real(readValue));
    EXPECT_FLOAT_EQ( std::imag(entry), std::imag(readValue));
}



}
