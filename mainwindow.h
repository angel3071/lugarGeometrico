#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <gsl/gsl_poly.h>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_calculate_clicked();

private:
    Ui::MainWindow *ui;
    bool itsCalculated = false;

};

#endif // MAINWINDOW_H
