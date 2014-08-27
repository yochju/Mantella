#ifndef SIDEWINDOW_HPP
#define SIDEWINDOW_HPP

#include <QDebug>

#include <QMainWindow>
#include <QVBoxLayout>
#include <QDialogButtonBox>
#include <QPushButton>

class DraggableTabWidget;

class SideWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit SideWindow(QMainWindow *parent = 0);
    ~SideWindow();

    QWidget *placeboWidget;
    QMenuBar *menuBar;
    DraggableTabWidget *tabWidget;

    DraggableTabWidget* getTabWidget();
signals:

public slots:

};

#endif // SIDEWINDOW_HPP
