#include <stdlib.h>
#include <stdio.h>
#include <iostream>
 
int int_sorter( const void *first_arg, const void *second_arg )
{
    int first = *(int*)first_arg;
    int second = *(int*)second_arg;
    if ( first < second )
    {
        return -1;
    }
    else if ( first == second )
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

class Sorter{
    public:
        virtual int compare (const void* first, const void* second)=0;
};

class AscendSorter: public Sorter {
    
    virtual int compare(const void* first_arg, const void* second_arg){
        int first = *(int*) first_arg;
        int second = *(int*) second_arg;

        if (first < second){
            return -1; 
        }
        else if (first == second){
            return 0;
        }
        else{
            return 1;
        }

    }

};

struct Test1{ 
    int _lol; 
    float _fun;
    Test1()=default;
    Test1(int lol, float fun) : _lol(lol), _fun(fun){};
    ~Test1()=default;
};
struct Test2{ 
    int _cool; 
    float _dry; 
    float _stop;
    Test2()=default;
    Test2(int cool, float dry, float stop) : _cool(cool), _dry(dry), _stop(stop){};
    ~Test2()=default;
};
struct Overall { 
    Test1* _test1; 
    Test2* _test2;
    Overall(Test1* test1){
        //Test1(test1->_lol,test1->_fun);
        this->_test1=test1;
        this->_test2=0;
    }
    Overall(Test2* test2){
        //Test2(test2->_cool, test2->_dry, test2->_stop);
        this->_test2=test2;
        this->_test1=0;
    }
};

void testFunct(Overall overall){
    if (overall._test1 != 0 && overall._test2 == 0){ 
        std::cout << "test1" << "\n";
        std::cout << overall._test1->_fun << "  " << overall._test1->_lol << "\n";
        }
    if (overall._test2 != 0 && overall._test1 == 0){ 
        std::cout << "test2" << "\n";
        std::cout << overall._test2->_cool << overall._test2->_dry << overall._test2->_stop << "\n";
        }
}

int main()
{
    //Test1 test1Obj(1,3.4);
    Test2 test2Obj(2,4.5,2.1);
    Overall overallObj(&test2Obj);
    testFunct(overallObj);

    int arr1[10];
    for ( int i = 0; i < 10; ++i )
    {
        arr1[ i ] = 10 - i;
    }
    void* ptrConv = static_cast<void*>(arr1);
    std::cout << arr1 << "  " << ptrConv << "\n";
    AscendSorter* AscendSorterPtr = new AscendSorter();
    Sorter* ptrSorter = AscendSorterPtr;
    printf( "%d\n", ptrSorter->compare(ptrConv,ptrConv) );

    int array[10];
    int i;
    /* fill array */
    for ( i = 0; i < 10; ++i )
    {
        array[ i ] = 10 - i;
    }
    qsort( array, 10 , sizeof( int ), &int_sorter );
    for ( i = 0; i < 10; ++i )
    {
        printf ( "%d\n" ,array[ i ] );
    }
    delete AscendSorterPtr;
 
}