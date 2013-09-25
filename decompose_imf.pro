#-------------------------------------------------
#
# Project created by QtCreator 2013-07-19T09:27:00
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QMAKE_CXXFLAGS += -std=c++0x -pedantic

TARGET = decompose_imf
TEMPLATE = app

PRECOMPILED_HEADER = stdafx.h
HEADERS  += stdafx.h \
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
    cpp_utils/visitor.h \
    qt_utils/serialize_props.h \
    cpp_utils/user_parameter.h \
    cpp_utils/cow_ptr.h \
    cpp_utils/bwt.h

SOURCES += main.cpp\
    gui_main_window.cpp \
    calculations.cpp \
    qt_utils/serialize_props.cpp

FORMS    += \
    gui_main_window.ui

LIBS += -L/usr/local/lib/ -lopencv_core -lopencv_imgproc -lopencv_highgui

OTHER_FILES +=
