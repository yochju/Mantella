#include "optionsTab.hpp"

OptionsTab::OptionsTab(QWidget *parent) :
    QWidget(parent)
{
    layout = new QGridLayout;

    onlineOptLabel = new QLabel("Online Optimisation");
    onlineOptCheckBox = new QCheckBox;
    evalLabel = new QLabel("Evaluation method");
    evalDropdown = new QComboBox;
    populateEval();

    experimentalSetupGroup = new QGroupBox("Experimental setup");
    optProblemLabel = new QLabel("Optimisation problem");
    optProblemDropdown = new QComboBox;
    optAlgoLabel = new QLabel("Optimisation algorithm");
    optAlgoDropdown = new QComboBox;

    generalOptionsGroup = new QGroupBox("General options");
    parallelLabel = new QLabel("Parallelisation");
    parallelCheckBox = new QCheckBox;

    specificOptions = constructSpecificOptions(NULL,NULL);

    //creating the buttonbar on the bottom and connecting the buttons
    startButton = new QPushButton("Start");
    startButton->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_R));
    stopAllButton = new QPushButton("Stop All");
    stopAllButton->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_S));
    buttonBox = new QDialogButtonBox();
    buttonBox->addButton(startButton,QDialogButtonBox::YesRole);
    buttonBox->addButton(stopAllButton,QDialogButtonBox::NoRole);
    connect(startButton,SIGNAL(clicked()),this,SLOT(startPressed()));
    connect(stopAllButton,SIGNAL(clicked()),this,SLOT(stopAllPressed()));

    //add all elements
    layout->addWidget(onlineOptLabel,0,0);
    layout->addWidget(onlineOptCheckBox,0,1);
    layout->addWidget(evalLabel,1,0);
    layout->addWidget(evalDropdown,1,1);

    //laying out setup group
    experimentalSetupLayout = new QVBoxLayout;
    experimentalSetupLayout->addWidget(optProblemLabel);
    experimentalSetupLayout->addWidget(optProblemDropdown);
    experimentalSetupLayout->addWidget(optAlgoLabel);
    experimentalSetupLayout->addWidget(optAlgoDropdown);
    experimentalSetupGroup->setLayout(experimentalSetupLayout);
    layout->addWidget(experimentalSetupGroup,2,0,5,2);

    //laying out general options group
    generalOptionsLayout = new QHBoxLayout;
    generalOptionsLayout->addWidget(parallelLabel);
    generalOptionsLayout->addWidget(parallelCheckBox);
    generalOptionsGroup->setLayout(generalOptionsLayout);
    layout->addWidget(generalOptionsGroup,7,0,2,2);

    layout->addWidget(specificOptions,9,0,1,2);

    layout->addWidget(buttonBox,11,0,1,2);

    setLayout(layout);
}

QWidget* OptionsTab::constructSpecificOptions(QString optProblemName, QString optAlgoName)
{
    QGroupBox *specOpsGroup = new QGroupBox("Specific options");
    QGridLayout *specOpsLayout = new QGridLayout;
    specOpsGroup->setLayout(specOpsLayout);
#ifdef QT_DEBUG
    QLabel *testy = new QLabel("testy");
    QComboBox *testBox = new QComboBox;
    QLabel *testy2 = new QLabel("testy2");
    QComboBox *testBox2 = new QComboBox;
    specOpsLayout->addWidget(testy);
    specOpsLayout->addWidget(testBox);
    specOpsLayout->addWidget(testy2);
    specOpsLayout->addWidget(testBox2);
#endif

    //if nothing is selected, return empty;
    if(optProblemName == NULL || optAlgoName == NULL || optProblemName.isEmpty() || optAlgoName.isEmpty()) {
        return specOpsGroup;
    }
}

void OptionsTab::populateEval()
{

}

void OptionsTab::populateOptProblem()
{

}

void OptionsTab::populateOptAlgo()
{

}

void OptionsTab::startPressed() {
    qDebug() << "start pressed";
    SideWindow *test = new SideWindow();
    test->show();
}

void OptionsTab::stopAllPressed() {
    qDebug() << "stopAll pressed";
}
