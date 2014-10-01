#include "draggabletabbar.hpp"
#include "draggabletabwidget.hpp"

DraggableTabBar::DraggableTabBar(QWidget *parent) :
    QTabBar(parent)
{
    setAcceptDrops(true);
    setSelectionBehaviorOnRemove (QTabBar::SelectLeftTab);
    setMovable (true);
}

void DraggableTabBar::dragEnterEvent(QDragEnterEvent *event)
{
    qDebug() << "dragenter event";
    const QMimeData *mimeData = event->mimeData();
    if(mimeData->formats().contains("action") && (mimeData->data("action") == "tab-dragging")) {
        event->accept();
    }
}

void DraggableTabBar::dragLeaveEvent(QDragLeaveEvent *event)
{
    qDebug() << "dragleave event";
    event->accept();
}

void DraggableTabBar::dragMoveEvent(QDragMoveEvent *event)
{
    qDebug() << "dragmove event";
    const QMimeData *mimeData = event->mimeData();
    if(mimeData->formats().contains("action") && (mimeData->data("action") == "tab-dragging")) {
        event->accept();
    }
}

void DraggableTabBar::dropEvent(QDropEvent *event)
{
    qDebug() << "drop event";

    const QMimeData *mimeData = event->mimeData();
    if(mimeData->formats().contains("action") && (mimeData->data("action") == "tab-dragging")) {
        MainWindow::dndHandled = true;
        //dont need to do anything if we are dropping on the same tabbar again
        if(event->source()==this) {
            return;
        }
        //dropping on a different tabbar!
        else {
            DraggableTabBar *sourceBar = (DraggableTabBar *)event->source();
            DraggableTabWidget *sourceTabWidget = (DraggableTabWidget *)sourceBar->parentWidget();
            QWidget *draggedTab = sourceTabWidget->widget(sourceBar->currentIndex());
            QIcon icon = sourceTabWidget->tabIcon(sourceBar->currentIndex());
            QString text = sourceTabWidget->tabText(sourceBar->currentIndex());
            DraggableTabWidget *targetTabWidget = (DraggableTabWidget *)parentWidget();
            sourceTabWidget->removeTab(sourceBar->currentIndex());
            targetTabWidget->addTab(draggedTab,icon,text);
        }
        event->accept();
    }
}

void DraggableTabBar::mouseMoveEvent(QMouseEvent *event)
{
    qDebug() << "mousemove event";
    /////making mouse enter and leave events cause qt doesnt do this when a mousebutton is pressed (wat)
    if(!lastPos.isNull()) {
        //entering
        if(geometry().contains(event->pos()) && !geometry().contains(lastPos)) {
            mouseEnterEvent(event);
        }
        //leaving
        if(!geometry().contains(event->pos()) && geometry().contains(lastPos)) {
            mouseLeaveEvent(event);
        }
    }
    lastPos=event->pos();
    /////end

    QTabBar::mouseMoveEvent(event);
}

void DraggableTabBar::mouseReleaseEvent(QMouseEvent *event)
{
    qDebug() << "mouserelease event";
    QTabBar::mouseReleaseEvent(event);
}

void DraggableTabBar::mousePressEvent(QMouseEvent *event)
{
    /////neeeded for mouse enter/leave
    lastPos = event->pos();
    /////end
    qDebug() << "mousepress event";

    QTabBar::mousePressEvent(event);
}

void DraggableTabBar::mouseLeaveEvent(QMouseEvent *event)
{
    qDebug() << "mouseleave event";

    //finishing the normal tabmoving
    QMouseEvent* finishMoveEvent = new QMouseEvent (QEvent::MouseMove, event->pos (), Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
    QTabBar::mouseReleaseEvent(finishMoveEvent);
    delete finishMoveEvent;
    finishMoveEvent = NULL;

    //starting the drag
    QDrag *drag = new QDrag(this);
    QMimeData *mimeData = new QMimeData;
    // a crude way to distinguish tab-dragging drops from other ones
    mimeData->setData("action", "tab-dragging");
    drag->setMimeData(mimeData);
    drag->exec();
    qDebug() << MainWindow::dndHandled;
    if(!MainWindow::dndHandled) {
        SideWindow *sideWindow = new SideWindow;
        sideWindow->show();
        DraggableTabBar *sourceBar = this;
        DraggableTabWidget *sourceTabWidget = (DraggableTabWidget *)sourceBar->parentWidget();
        QWidget *draggedTab = sourceTabWidget->widget(sourceBar->currentIndex());
        QIcon icon = sourceTabWidget->tabIcon(sourceBar->currentIndex());
        QString text = sourceTabWidget->tabText(sourceBar->currentIndex());
        DraggableTabWidget *targetTabWidget = (DraggableTabWidget *)sideWindow->getTabWidget();
        sourceTabWidget->removeTab(sourceBar->currentIndex());
        targetTabWidget->addTab(draggedTab,icon,text);
        QTabBar *tabBar = targetTabWidget->findChild<QTabBar *>();
        tabBar->setTabButton(0,QTabBar::RightSide,0);
        tabBar->setTabButton(0,QTabBar::LeftSide,0);
    }
    MainWindow::dndHandled = false;
}

void DraggableTabBar::mouseEnterEvent(QMouseEvent *event)
{
    qDebug() << "mouseenter event";
}
