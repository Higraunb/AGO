#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "Algorithm.h"
#include <qwt_plot.h>
#include <qwt_plot_grid.h>
#include <qwt_legend.h>
#include <qwt_plot_curve.h>
#include <qwt_symbol.h>
#include <qwt_plot_magnifier.h>
#include <qwt_plot_panner.h>
#include <qwt_plot_picker.h>
#include <qwt_picker_machine.h>
#include <qwt_plot_marker.h>  // Добавляем для маркеров осей

#include <functional>

struct FuncInfo
{
    double lowerBound;
    double upperBound;
    std::function<double(double)> function;
    QColor color;
    QString funcName;

    FuncInfo():lowerBound(0), upperBound(1),
        function([](double x) { return x; }),
        color(Qt::black), funcName("y = x")
    {}

    FuncInfo(const double &lowerBound_, const double &upperBound_,
             const std::function<double (double)> &function_,
             QColor color_, const QString &funcName_):
        lowerBound(lowerBound_),upperBound(upperBound_),
        function(function_), color(color_), funcName(funcName_)
    {
    }
};

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    FuncInfo function;
    QwtPlotCurve *curve;
    QPolygonF points;
    QPointF a;
    TAlgorithm<double,1> algo;
    QList<QwtPlotMarker*> pointMarkers;
    void addPoint(vector<double> other, QColor color = Qt::red, int size = 8);
    void addPlot();
    void addPlotGrid();
    void addCurve(const FuncInfo& func);
    void enableMagnifier();
    void enableMovingOnPlot();
    void enablePicker();
    void addAxes();
};

#endif // MAINWINDOW_H
