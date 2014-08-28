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
    setupMenuBar();
}


SideWindow::~SideWindow()
{
}

DraggableTabWidget* SideWindow::getTabWidget()
{
    return tabWidget;
}

void SideWindow::setupMenuBar()
{
    menuBar = new QMenuBar;
    setMenuBar(menuBar);
    //start menu
    QMenu *startMenu = menuBar->addMenu("Start");
    exitAct = new QAction("Exit",this);
    exitAct->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_Q));
    connect(exitAct,SIGNAL(triggered()),this,SLOT(exitProgram()));
    startMenu->addAction(exitAct);
    //view menu
    QMenu *viewMenu = menuBar->addMenu("View");
    //help menu
    QMenu *helpMenu = menuBar->addMenu("Help");
    aboutAct = new QAction("About",this);
    connect(aboutAct,SIGNAL(triggered()),this,SLOT(about()));
    helpMenu->addAction(aboutAct);

    //adding to menubar
    menuBar->addMenu(startMenu);
    menuBar->addMenu(viewMenu);
    menuBar->addMenu(helpMenu);
}

void SideWindow::about()
{
    QMessageBox::about(this,"im titling in your title","This is online optimisation. Yay!");
}

void SideWindow::exitProgram()
{
    QMessageBox msgBox;
    msgBox.setText("Are you sure you want to quit?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);
    int ret = msgBox.exec();
    if(ret==QMessageBox::Yes) {
        QCoreApplication::quit();
    }
}
