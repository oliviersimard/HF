//#include "../src/Hubbard.cpp"
#include <gtest/gtest.h>


// Integral1D tests
//static Integral1D int1D(2.0,0.1+0.2*im);
//static Integral1D int1D2(3.0,-0.1+0.2*im);
//// Integral2D tests
//static Integral2D int2D(2.0,1.0,0.1+0.2*im);
//static Integral2D int2D2(3.0,-1.0,-0.1+0.2*im);



TEST(Integral1DTest, ContainerTest1D){
    ASSERT_EQ(0, 0);
//    ASSERT_EQ(2.0,int1D.qx);
//    ASSERT_EQ(0.1+0.2*im,int1D.iwn);
//    ASSERT_EQ(5.0,(int1D+int1D2).qx);
//    ASSERT_EQ(0.0+0.4*im,(int1D+int1D2).iwn);
//    ASSERT_EQ(-1.0,(int1D-int1D2).qx);
//    ASSERT_EQ(0.2+0.0*im,(int1D-int1D2).iwn);
}

TEST(Integral2DTest, ContainerTest2D){
//    ASSERT_EQ(2.0,int2D.qx);
//    ASSERT_EQ(1.0,int2D.qy);
//    ASSERT_EQ(0.1+0.2*im,int2D.iwn);
//    ASSERT_EQ(5.0,(int2D+int2D2).qx);
//    ASSERT_EQ(0.0,(int2D+int2D2).qy);
//    ASSERT_EQ(0.0+0.4*im,(int2D+int2D2).iwn);
//    ASSERT_EQ(-1.0,(int2D-int2D2).qx);
//    ASSERT_EQ(2.0,(int2D-int2D2).qy);
//    ASSERT_EQ(0.2+0.0*im,(int2D-int2D2).iwn);
}


int main(int argc, char **argv) {

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
