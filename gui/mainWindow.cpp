#include "mainWindow.hpp"
#include "draggabletabwidget.hpp"

bool MainWindow::dndHandled = false;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{


    tabWidget = new DraggableTabWidget;
    //necessary to prevent optionstab from being closed
    tabWidget->addTab(new OptionsTab,"Start");
    QTabBar *tabBar = tabWidget->findChild<QTabBar *>();
    tabBar->setTabButton(0,QTabBar::RightSide,0);
    tabBar->setTabButton(0,QTabBar::LeftSide,0);

    //gotta have a layout, the placebowidget is needed for weird qt reasons (cant set layout of mainwindow)
    placeboWidget = new QWidget;
    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(tabWidget);
    placeboWidget->setLayout(mainLayout);
    setCentralWidget(placeboWidget);

    setWindowTitle("Online Optimisation");
    setupMenuBar();

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

void MainWindow::setupMenuBar()
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

void MainWindow::about()
{
    QMessageBox::about(this,"im titling in your title","This is online optimisation. Yay!");
}

void MainWindow::exitProgram()
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

