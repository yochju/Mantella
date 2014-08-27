#ifndef OPTIONSTAB_HPP
#define OPTIONSTAB_HPP

#include <QDebug>

#include <QWidget>
#include <QLabel>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLineEdit>
#include <QGroupBox>
#include <QComboBox>
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QPushButton>
#include "sidewindow.hpp"

class OptionsTab : public QWidget
{
    Q_OBJECT
public:
    explicit OptionsTab(QWidget *parent = 0);

    QGridLayout *layout;

    QLabel *onlineOptLabel;
    QCheckBox *onlineOptCheckBox;
    QLabel *evalLabel;
    QComboBox *evalDropdown;
    QGroupBox *experimentalSetupGroup;
    QVBoxLayout *experimentalSetupLayout;
    QLabel *optProblemLabel;
    QComboBox *optProblemDropdown;
    QLabel *optAlgoLabel;
    QComboBox *optAlgoDropdown;
    QGroupBox *generalOptionsGroup;
    QHBoxLayout *generalOptionsLayout;
    QLabel *parallelLabel;
    QCheckBox *parallelCheckBox;
    QWidget *specificOptions;

    QWidget* constructSpecificOptions(QString optProblemName, QString optAlgoName);
    void populateEval();
    void populateOptProblem();
    void populateOptAlgo();

    QDialogButtonBox *buttonBox;
    QPushButton *startButton;
    QPushButton *stopAllButton;

signals:

public slots:
    void startPressed();
    void stopAllPressed();

};

#endif // OPTIONSTAB_H
