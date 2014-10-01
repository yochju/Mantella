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
    populateOptAlgo();

    generalOptionsGroup = new QGroupBox("General options");
    parallelLabel = new QLabel("Parallelisation");
    parallelCheckBox = new QCheckBox;
    populateOptProblem();

    specificOptions = constructSpecificOptions(OptimisationProblem::Type::AttractiveSectorFunction);

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

    layout->addWidget(buttonBox,10,0,1,2);

    setLayout(layout);
}

QWidget* OptionsTab::constructSpecificOptions(OptimisationProblem::Type optProblem)
{
    //create a "container"
    QGroupBox *specOpsGroup = new QGroupBox("Specific options");
    QGridLayout *specOpsLayout = new QGridLayout;
    specOpsGroup->setLayout(specOpsLayout);

    //DEBUG STUFF
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

    //fill "container" depending on problem
    switch(optProblem) {
    case OptimisationProblem::Type::AttractiveSectorFunction:
    {
        QLabel *nametest = new QLabel(QString::fromStdString(OptimisationProblem::getName(optProblem)));
        specOpsLayout->addWidget(nametest);
        break;
    }
    case OptimisationProblem::Type::BentCigarFunction:
    {
        QLabel *nametest = new QLabel(QString::fromStdString(OptimisationProblem::getName(optProblem)));
        specOpsLayout->addWidget(nametest);
        QLabel *nametest2 = new QLabel(QString::fromStdString(OptimisationProblem::getName(optProblem)));
        specOpsLayout->addWidget(nametest2);
        QLabel *nametest3 = new QLabel(QString::fromStdString(OptimisationProblem::getName(optProblem)));
        specOpsLayout->addWidget(nametest3);
        QLabel *nametest4 = new QLabel(QString::fromStdString(OptimisationProblem::getName(optProblem)));
        specOpsLayout->addWidget(nametest4);
        QLabel *nametest5 = new QLabel(QString::fromStdString(OptimisationProblem::getName(optProblem)));
        specOpsLayout->addWidget(nametest5);
        break;
    }
    case OptimisationProblem::Type::BuecheRastriginFunction:

    case OptimisationProblem::Type::CompositeGriewankRosenbrockFunctionF8F2:

    case OptimisationProblem::Type::DifferentPowersFunction:

    case OptimisationProblem::Type::DiscusFunction:

    case OptimisationProblem::Type::EllipsoidalFunction:

    case OptimisationProblem::Type::EllipsoidalFunctionRotated:

    case OptimisationProblem::Type::GallaghersGaussian101mePeaksFunction:

    case OptimisationProblem::Type::GallaghersGaussian21hiPeaksFunction:

    case OptimisationProblem::Type::KatsuuraFunction:

    case OptimisationProblem::Type::LinearSlope:

    case OptimisationProblem::Type::LunacekBiRastriginFunction:

    case OptimisationProblem::Type::RastriginFunction:

    case OptimisationProblem::Type::RastriginFunctionRotated:

    case OptimisationProblem::Type::RosenbrockFunction:

    case OptimisationProblem::Type::RosenbrockFunctionRotated:

    case OptimisationProblem::Type::SchaffersF7Function:

    case OptimisationProblem::Type::SchaffersF7FunctionIllConditioned:

    case OptimisationProblem::Type::SchwefelFunction:

    case OptimisationProblem::Type::SharpRidgeFunction:

    case OptimisationProblem::Type::SphereFunction:

    case OptimisationProblem::Type::StepEllipsoidalFunction:

    case OptimisationProblem::Type::WeierstrassFunction:

        break;
    }

    return specOpsGroup;
}

void OptionsTab::populateEval()
{

}

Q_DECLARE_METATYPE(OptimisationProblem::Type)

void OptionsTab::populateOptProblem()
{
    OptimisationProblem::_names.at(OptimisationProblem::Type::AttractiveSectorFunction);
    for(auto &any : OptimisationProblem::_names) {
        qDebug() << typeid(any.first).name();
        optProblemDropdown->addItem(QString::fromStdString(any.second),QVariant::fromValue(any.first));
    }
    connect(optProblemDropdown,SIGNAL(currentIndexChanged(int)),this,SLOT(optProblemSelectionChanged(int)));
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

void OptionsTab::optAlgoSelectionChanged(int)
{

}

void OptionsTab::optProblemSelectionChanged(int)
{
    QVariant var = optProblemDropdown->currentData();
    OptimisationProblem::Type problemType = var.value<OptimisationProblem::Type>();

    QWidget *updatedSpecificOptions = constructSpecificOptions(problemType);
    layout->removeWidget(specificOptions);
    specificOptions->deleteLater();
    layout->addWidget(updatedSpecificOptions,9,0,1,2);
    specificOptions = updatedSpecificOptions;

}

void OptionsTab::evalSelectionChanged(int)
{

}
