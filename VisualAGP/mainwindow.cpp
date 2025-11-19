#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <cmath>

double Func(double x)
{
    return sin(x);
}

double Func1(const TPoint<double, 1>& other)
{
    double x = other[0];
    return sin(x);
}
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    function(),
    curve(nullptr),
    algo(),
    points(),
    a(0, 0)
{
    ui->setupUi(this);
    algo = TAlgorithm<double , 1>(TPoint<double, 1>(-2.0),TPoint<double, 1>(2.0), 0.01, 2, Func1);
    algo.AGPStronginaMin();
    function = FuncInfo(-2.0, 5.0,
                        Func,
                        Qt::blue,
                        "y = 2x - 1");
    ui->textEditInform->setReadOnly(true);
    ui->textEditInform->setFontPointSize(14);
    ui->textEditInform->setText(QString::fromStdString(algo.GetInfo()));
    addPlot();
    addPlotGrid();
    addAxes();
    addCurve(function);
    addPoint(algo.GetAllPoint());
    enableMagnifier();
    enableMovingOnPlot();
    enablePicker();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::addPlot()
{
    ui->PlotWidget->setTitle("Function Plot");
    ui->PlotWidget->setAxisTitle(QwtPlot::xBottom, "x");
    ui->PlotWidget->setAxisTitle(QwtPlot::yLeft, "y");
}

void MainWindow::addPlotGrid()
{
    if (!ui->PlotWidget) return;

    QwtPlotGrid *grid = new QwtPlotGrid();
    grid->setPen(Qt::gray, 0.0, Qt::DotLine);
    grid->attach(ui->PlotWidget);
}

void MainWindow::addAxes()
{
    if (!ui->PlotWidget) return;

    // Создаем маркер для оси Y (x = 0)
    QwtPlotMarker *yAxisMarker = new QwtPlotMarker();
    yAxisMarker->setLineStyle(QwtPlotMarker::VLine);  // Вертикальная линия
    yAxisMarker->setLinePen(QPen(Qt::black, 2, Qt::SolidLine));  // Черная, толщина 2
    yAxisMarker->setXValue(0.0);  // Позиция x = 0
    yAxisMarker->attach(ui->PlotWidget);

    // Создаем маркер для оси X (y = 0)
    QwtPlotMarker *xAxisMarker = new QwtPlotMarker();
    xAxisMarker->setLineStyle(QwtPlotMarker::HLine);  // Горизонтальная линия
    xAxisMarker->setLinePen(QPen(Qt::black, 2, Qt::SolidLine));  // Черная, толщина 2
    xAxisMarker->setYValue(0.0);  // Позиция y = 0
    xAxisMarker->attach(ui->PlotWidget);
}

void MainWindow::addPoint(vector<double> other, QColor color, int size)
{
    if (!ui->PlotWidget) return;
    double x = 0;
    double y = 0;
    for (auto& point : other)
    {
        x = point;
        // Создаем маркер для точки
        QwtPlotMarker *pointMarker = new QwtPlotMarker();

        // Устанавливаем позицию точки
        pointMarker->setValue(x, y);

        // Создаем символ (круг)
        QwtSymbol *symbol = new QwtSymbol(QwtSymbol::Ellipse);
        symbol->setSize(size, size);           // Размер
        symbol->setColor(color);               // Цвет заливки
        symbol->setPen(QPen(Qt::black, 1));    // Черная обводка

        // Устанавливаем символ для маркера
        pointMarker->setSymbol(symbol);

        // Прикрепляем маркер к графику
        pointMarker->attach(ui->PlotWidget);
        pointMarkers.append(pointMarker);
    }
    ui->PlotWidget->replot();
}

void MainWindow::addCurve(const FuncInfo& func)
{
    if (!ui->PlotWidget) return;

    // Очищаем предыдущие точки
    points.clear();

    // Генерируем точки для функции
    int pointCount = 200;
    for (int i = 0; i <= pointCount; ++i) {
        double x = func.lowerBound + (func.upperBound - func.lowerBound) * i / pointCount;
        double y = func.function(x);
        points << QPointF(x, y);
    }

    // Создаем кривую
    curve = new QwtPlotCurve(func.funcName);
    curve->setSamples(points);
    curve->setPen(func.color, 2);
    curve->setRenderHint(QwtPlotItem::RenderAntialiased);
    curve->attach(ui->PlotWidget);

    // Перерисовываем график
    ui->PlotWidget->replot();
}

void MainWindow::enableMagnifier()
{
    if (!ui->PlotWidget) return;

    QwtPlotMagnifier *magnifier = new QwtPlotMagnifier(ui->PlotWidget->canvas());
    magnifier->setMouseButton(Qt::MidButton);
}

void MainWindow::enableMovingOnPlot()
{
    if (!ui->PlotWidget) return;

    QwtPlotPanner *panner = new QwtPlotPanner(ui->PlotWidget->canvas());
    panner->setMouseButton(Qt::LeftButton);
}

void MainWindow::enablePicker()
{
    if (!ui->PlotWidget) return;

    QwtPlotPicker *picker = new QwtPlotPicker(QwtPlot::xBottom, QwtPlot::yLeft,
                                              QwtPicker::CrossRubberBand,
                                              QwtPicker::AlwaysOn,
                                              ui->PlotWidget->canvas());
    picker->setStateMachine(new QwtPickerDragPointMachine());
}
