#pragma once

#include <QDebug>

#include <QTabWidget>
#include <QtGui>
#include <QShortcut>
#include <QAction>
#include <QMessageBox>
#include "draggabletabbar.hpp"

class DraggableTabWidget : public QTabWidget
{
    Q_OBJECT
public:
    explicit DraggableTabWidget(QWidget *parent = 0);

    virtual void removeTab(int index);

signals:

public slots:
    //shortcuts
    void closeTab(int index = -2);

protected:
    virtual void dragEnterEvent(QDragEnterEvent *event);
    virtual void dragLeaveEvent(QDragLeaveEvent *event);
    virtual void dragMoveEvent(QDragMoveEvent *event);
    virtual void dropEvent(QDropEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void mousePressEvent(QMouseEvent *event);

};
