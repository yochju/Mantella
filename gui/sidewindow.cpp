#include "sidewindow.hpp"
#include "draggabletabwidget.hpp"

SideWindow::SideWindow(QMainWindow *parent) :
    QMainWindow(parent)
{
    tabWidget = new DraggableTabWidget;

    //gotta have a layout, the placebowidget is needed for weird qt reasons (cant set layout of mainwindow)
    placeboWidget = new QWidget;
    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(tabWidget);
    placeboWidget->setLayout(mainLayout);
    setCentralWidget(placeboWidget);

    setWindowTitle("Online Optimisation");
}


SideWindow::~SideWindow()
{
}

DraggableTabWidget* SideWindow::getTabWidget()
{
    return tabWidget;
}
