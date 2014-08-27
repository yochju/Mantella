#include "mainWindow.hpp"
#include "draggabletabwidget.hpp"

bool MainWindow::dndHandled = false;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{


    tabWidget = new DraggableTabWidget;
    //necessary to prevent optionstab from being closed
    tabWidget->addTab(new OptionsTab,"Start");
//    QTabBar *tabBar = tabWidget->findChild<QTabBar *>();
//    tabBar->setTabButton(0,QTabBar::RightSide,0);
//    tabBar->setTabButton(0,QTabBar::LeftSide,0);

    //gotta have a layout, the placebowidget is needed for weird qt reasons (cant set layout of mainwindow)
    placeboWidget = new QWidget;
    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(tabWidget);
    placeboWidget->setLayout(mainLayout);
    setCentralWidget(placeboWidget);

    setWindowTitle("Online Optimisation");
#ifdef QT_DEBUG
    addTestTabs();
#endif
}

MainWindow::~MainWindow()
{
}

void MainWindow::addTestTabs() {
    for(int i = 0; i < 5; i++) {
        tabWidget->addTab(new QWidget,QString("lol")+QString::number(i));
    }
}

void MainWindow::startPressed() {
    qDebug() << "start pressed";
    SideWindow *test = new SideWindow();
    test->show();
}

void MainWindow::stopAllPressed() {
    qDebug() << "stopAll pressed";
}

