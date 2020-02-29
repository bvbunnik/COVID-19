TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -L/usr/local/lib -lgsl -lgslcblas -lm

SOURCES += \
        main.cpp \
    solvebeta.cpp

DISTFILES += \
    parameters_used.txt

HEADERS += \
    solvebeta.h
