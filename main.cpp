#include "gui_main_window.h"
#include "qt_utils/exception_handling_application.h"

int main(int argc, char *argv[])
{
    qu::ExceptionHandlingApplication a(argc, argv);
    gui::MainWindow w;
    w.show();
    
    return a.exec();
}
