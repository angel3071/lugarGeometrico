#-------------------------------------------------
#
# Project created by QtCreator 2015-11-03T23:42:43
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Pruebita
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    qcustomplot.cpp

HEADERS  += mainwindow.h \
    qcustomplot.h

FORMS    += mainwindow.ui

DISTFILES +=
LIBS += -L/usr/local/lib -lgsl -lgslcblas -lm
