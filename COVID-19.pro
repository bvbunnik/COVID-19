TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lgsl -lgslcblas -lm

SOURCES += \
        main.cpp

DISTFILES += \
    parameters_used.txt
