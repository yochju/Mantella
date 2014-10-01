#pragma once

#include <QDebug>

#include <QVBoxLayout>
#include <QMainWindow>
#include <QTabBar>
#include <QMenuBar>
#include <QMenu>
#include <QAction>
#include <QMessageBox>
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


    //used by the menu
    void setupMenuBar();
    QAction *exitAct;
    QAction *aboutAct;

public slots:
    //used by the menu
    void about();
    void exitProgram();
};

