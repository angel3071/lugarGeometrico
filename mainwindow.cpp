#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <gsl/gsl_poly.h>
#include <QPen>
#include <QtMath>

// Este arreglo de plumas permite que se vayan intercalando los colores de las trazas
QPen pens[10] = {QPen(Qt::blue),QPen(Qt::red),QPen(Qt::magenta),QPen(Qt::green),QPen(Qt::black),
                 QPen(Qt::gray),QPen(Qt::yellow),QPen(Qt::darkBlue),QPen(Qt::darkRed),QPen(Qt::darkYellow)};

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_calculate_clicked()
{
    if(itsCalculated)return;
    int gradeDS = 5;
    int gradeNS = 1;
    double k = 0;
    //double ns[gradeNS + 1] = {1};
    double ns[gradeNS + 1] = {3,1}; // k(s+3)
    double gs[gradeDS + 1] = {0,60,82,54,13,1}; // G(s) = s^5 + 13s^4 + 54s^3 + 82s^2 + 60s + 0
    //double gs[gradeDS + 1] = {0,2,3,1};
    double newgs[gradeDS + 1] = {0,60,82,54,13,1};
    double zns[gradeNS * 2];
/// Raices en malla abierta
    double zgs[gradeDS * 2];
    bool isComplex[gradeDS + 1];
    //Se prepara el espacio en memoria para la solución del polinomio y se resuelve
    gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(gradeDS + 1);
    gsl_poly_complex_solve (gs, gradeDS + 1, w, zgs);
    gsl_poly_complex_workspace_free (w);
    w = gsl_poly_complex_workspace_alloc(gradeNS + 1);
    if(gradeNS!=0){
        gsl_poly_complex_solve (ns, gradeNS + 1, w, zns);
        gsl_poly_complex_workspace_free (w); //Se libera el espacio en memoria
    }
/// Polos en el plano, se agregan a la grafica
    QCustomPlot *customPlot = ui->customplot; //El objeto donde se grafica
    for(int i=0; i<gradeDS; i++){
        customPlot->addGraph();
        customPlot->graph(i)->setData(QVector<double>() << zgs[2*i], QVector<double>() << zgs[2*i+1]);
        customPlot->graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCross, 5));
        customPlot->graph(i)->setName(tr("p").append(QString::number(i)));
    }
/// Zeros en el plano, se agregan a la gráfica
    int graphs = customPlot->graphCount();
    int d = 0;
    for(int i = graphs; i < graphs + gradeNS; i++){
        customPlot->addGraph();
        customPlot->graph(i)->setData(QVector<double>() << zns[2*d], QVector<double>() << zns[2*d+1]);
        customPlot->graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
        customPlot->graph(i)->setName(tr("z").append(QString::number(d)));
        d++;
    }
/// Polos, zeros, centroide y grado relativo en panel, solo de manera informativa
    this->ui->results->setText(tr("Polos:\n\n"));
    for(int i=0; i<gradeDS; i++)
        this->ui->results->setText(ui->results->toPlainText().append("p").append(QString::number(i))
                                   .append(" = ").append(QString::number( zgs[2*i],'g',3).append(" + ")
                                   .append(QString::number(zgs[2*i+1],'g',3))).append("i").append("\n\n"));
    this->ui->results->setText(ui->results->toPlainText().append(tr("Zeros:\n\n")));
    for(int i=0; i<gradeNS; i++)
        this->ui->results->setText(ui->results->toPlainText().append("z").append(QString::number(i))
                                   .append(" = ").append(QString::number( zns[2*i],'g',3).append(" + ")
                                   .append(QString::number(zns[2*i+1],'g',3))).append("i").append("\n\n"));
    int GR = gradeDS - gradeNS;
    this->ui->results->setText(this->ui->results->toPlainText().append(tr("Grado Relativo:  ")));
    this->ui->results->setText(this->ui->results->toPlainText().append(QString::number(GR)));

    double centroide = 0.0;
    for(int i=0;i<gradeDS;i++){
        centroide+=zgs[2*i];
    }
    double tmp = 0.0;
    for(int i=0;i<gradeNS;i++)
        tmp+=zns[2*i];
    centroide-=tmp;centroide/=GR;
    this->ui->results->setText(this->ui->results->toPlainText().append(tr("\n\n")));
    this->ui->results->setText(this->ui->results->toPlainText().append(tr("Centroide:  ")));
    this->ui->results->setText(this->ui->results->toPlainText().append(QString::number(centroide)));

    QVector<double> rootsReal[gradeDS], rootsImg[gradeDS];
    for(int i =0; i< gradeDS; i++){
        rootsReal[i] = QVector<double>(4000);
        rootsImg[i] = QVector<double>(4000);
    }

/// Traza del lugar geométrico
    for(k = 0; k < 4000; k++){

        for(int i =0; i<gradeNS + 1; i++){
            newgs[i] = gs[i] + ((double)k-1999/(10.0))*ns[i]; //Incrementos de 0.1 en k, desde -2000 hasta 2000
        }
        w = gsl_poly_complex_workspace_alloc(gradeDS + 1);
        gsl_poly_complex_solve(newgs, gradeDS + 1, w, zgs);
        gsl_poly_complex_workspace_free (w);
        for(int i =0;i<gradeDS;i++){
            if(zgs[2*i+1] != 0.0)
                isComplex[i] = true;

            if(isComplex[i] && zgs[2*i+1] != 0.0){
                rootsReal[i][k] = zgs[2*i];
                rootsImg[i][k] = zgs[2*i+1];
            }
            qDebug() << "Re " << QString::number(rootsReal[i][k]) << " Im "
                     << QString::number(rootsImg[i][k]) << " k " << QString::number(((double)k/(10.0))) << " pol "
                     << QString::number(newgs[0]) << ","<< QString::number(newgs[1]) << ","<< QString::number(newgs[2]) << ","
                     << QString::number(newgs[3]) << ","<< QString::number(newgs[4]) << ","<< QString::number(newgs[5]);

        }


    }

    qDebug() << "graph count " + QString::number(customPlot->graphCount());

    d =0;
    graphs = customPlot->graphCount();
    for(int i =graphs;i<graphs + gradeDS; i ++){
        customPlot->addGraph();
        customPlot->graph(i)->setData(rootsReal[d], rootsImg[d]);
        customPlot->graph(i)->setName(tr("z").append(QString::number(d)));
        customPlot->graph(i)->setPen(pens[d]);
        d++;
    }
    d=0;
    if(GR>0){
        double angle;
        graphs = customPlot->graphCount();
        for(int i =graphs; i < GR + graphs; i++){
            angle = 0.0174533*(180.0*(double)(2*d + 1))/(double)GR;
            angle = qTan(angle);
            d++;
        }
    }

/// Configuraciones del plot
    qDebug() << "graph count " + QString::number(customPlot->graphCount());
    // give the axes some labels:
    customPlot->xAxis->setLabel("Re");
    customPlot->yAxis->setLabel("Im");
    // set axes ranges, so we see all data:
    customPlot->rescaleAxes();
    customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    customPlot->legend->setVisible(true);
    customPlot->legend->setBrush(QBrush(QColor(255,255,255,150)));
    customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);
    customPlot->replot();



    itsCalculated = true;

}
