#-------------------------------------------------
#
# Project created by QtCreator 2018-04-01T17:19:31
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = VehicleControlPredict_Symmetrical
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
    mainwindow.cpp \
    main.cpp \
    integrate.cpp \
    qcustomplot.cpp \
    predict.cpp \
    symm_functions.cpp \
    DOPRI8_symmetrical_plot.cpp \
    DOPRI8_final_plot.cpp \
    DOPRI8_custom.cpp \
    custom_functions.cpp \
    control_solve.cpp \
    custom_control_find.cpp

HEADERS += \
    mainwindow.h \
    include.h \
    qcustomplot.h \
    vector.h \
    custom_functions.h \
    symm_functions.h \
    parameters.h

FORMS += \
        mainwindow.ui

DISTFILES +=

SUBDIRS += \
    VehicleControlPredict_Symmetrical.pro

LIBS += -L/usr/include/gsl -lgsl -lgslcblas -lm
