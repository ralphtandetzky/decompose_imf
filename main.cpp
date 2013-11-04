/*#include "gui_main_window.h"
#include "qt_utils/exception_handling_application.h"


int main(int argc, char *argv[])
{
    qu::ExceptionHandlingApplication a(argc, argv);
    gui::MainWindow w;
    w.show();
    
    return a.exec();
}
*/

#include "cpp_utils/cow_tree.h"
#include <iostream>

int main()
try
{
    cu::CowMap<int,int> t;

    static const int N = 100;

    for ( int i = 0; i < N; ++i )
        t.insert( i, N-i );

    for ( int i = 0; i < N; ++i )
        t.modify( i, []( int & i ){ i = i*i; } );

    for ( int i = N/2; i < N; ++i )
        t.remove( i );

    t.modifyAll( [](int,int&i){++i;} );

    t.readAll( [](int key, int val)
    {
        std::cout
            << "key = "    << key
            << ",\tval = " << val
            << std::endl;
    });
}
catch( std::exception & e )
{
    std::cout << e.what() << std::endl;
}
