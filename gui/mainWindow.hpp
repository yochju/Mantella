#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QDebug>

#include <QVBoxLayout>
#include <QMainWindow>
#include <QTabBar>
#include "optionsTab.hpp"
#include "sidewindow.hpp"

class DraggableTabWidget;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    QWidget *placeboWidget;
    QMenuBar *menuBar;
    DraggableTabWidget *tabWidget;

    void addTestTabs();

    static bool dndHandled;
    static bool dndInside;
    static QWidget draggedObject;

public slots:
    void startPressed();
    void stopAllPressed();

};

#endif // MAINWINDOW_HPP
