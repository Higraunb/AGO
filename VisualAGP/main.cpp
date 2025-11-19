#include <QApplication>
#include "mainwindow.h"
#include <QDebug>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    qDebug() << "Starting application with Qwt test...";

    try {
        MainWindow w;
        w.show();
        w.setWindowTitle("Qwt Library Test - VisualAGP");

        qDebug() << "Qwt test application started successfully!";

        return a.exec();

    } catch (const std::exception& e) {
        qCritical() << "Failed to start application:" << e.what();
        return -1;
    }
}
