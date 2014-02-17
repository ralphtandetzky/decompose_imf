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
    cpp_utils/bwt.h \
    qt_utils/gui_user_parameter.h \
    cpp_utils/user_parameter_container.h \
    cpp_utils/cloning.h \
    qt_utils/gui_property_sheet.h \
    cpp_utils/formula_parser.h \
    cpp_utils/exception.h \
    cpp_utils/exception_handling.h \
    qt_utils/exception_handling.h \
    qt_utils/exception_handling_application.h \
    cpp_utils/cow_map.h \
    qt_utils/event_filter.h \
    qt_utils/invoke_in_thread.h \
    cpp_utils/extract_by_line.h \
    processing.h

SOURCES += main.cpp\
    gui_main_window.cpp \
    calculations.cpp \
    qt_utils/serialize_props.cpp \
    qt_utils/gui_user_parameter.cpp \
    cpp_utils/user_parameter_container.cpp \
    cpp_utils/user_parameter.cpp \
    qt_utils/gui_property_sheet.cpp \
    cpp_utils/formula_parser.cpp \
    qt_utils/exception_handling.cpp \
    cpp_utils/extract_by_line.cpp \
    processing.cpp

FORMS    += \
    gui_main_window.ui

LIBS += -L/usr/local/lib/ -lopencv_core -lopencv_imgproc -lopencv_highgui

OTHER_FILES += \
    notes.txt

RESOURCES += \
    resources.qrc
