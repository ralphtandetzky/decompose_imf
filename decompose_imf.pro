#-------------------------------------------------
#
# Project created by QtCreator 2013-07-19T09:27:00
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QMAKE_CXXFLAGS += -std=c++0x

TARGET = decompose_imf
TEMPLATE = app

SOURCES += main.cpp\
    gui_main_window.cpp \
    calculations.cpp

HEADERS  += \
    calculations.h \
    gui_main_window.h \
    cpp_utils/concurrent_queue.h \
    cpp_utils/locking.h \
    cpp_utils/math_constants.h \
    cpp_utils/more_algorithms.h \
    cpp_utils/optimize.h \
    cpp_utils/parallel_executor.h \
    cpp_utils/scope_guard.h \
    cpp_utils/spin_lock.h \
    cpp_utils/sqr.h \
    cpp_utils/std_make_unique.h \
    cpp_utils/value_ptr.h \
    cpp_utils/virtual_call.h \
    cpp_utils/visitor.h

FORMS    += \
    gui_main_window.ui

LIBS += -L/usr/local/lib/ -lopencv_core -lopencv_imgproc -lopencv_highgui

OTHER_FILES +=
