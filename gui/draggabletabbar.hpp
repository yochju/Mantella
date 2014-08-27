#ifndef DRAGGABLETABBAR_HPP
#define DRAGGABLETABBAR_HPP

#include <QDebug>

#include <QTabBar>
#include <QtGui>
#include "mainWindow.hpp"
#include "sidewindow.hpp"

class DraggableTabWidget;

class DraggableTabBar : public QTabBar
{
    Q_OBJECT
public:
    explicit DraggableTabBar(QWidget *parent = 0);

signals:

public slots:


protected:
    virtual void dragEnterEvent(QDragEnterEvent *event);
    virtual void dragLeaveEvent(QDragLeaveEvent *event);
    virtual void dragMoveEvent(QDragMoveEvent *event);
    virtual void dropEvent(QDropEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void mousePressEvent(QMouseEvent *event);

    /////needed for mouse enter and leave
    virtual void mouseLeaveEvent(QMouseEvent *event);
    virtual void mouseEnterEvent(QMouseEvent *event);
    QPoint lastPos;
};

#endif // DRAGGABLETABBAR_HPP
