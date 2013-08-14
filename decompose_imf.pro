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
    gui_main_window.h \
    optimize.h \
    calculations.h

FORMS    += \
    gui_main_window.ui

LIBS += -L/usr/local/lib/ -lopencv_core -lopencv_imgproc -lopencv_highgui
