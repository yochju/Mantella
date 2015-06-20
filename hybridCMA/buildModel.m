function model = buildModel(cur_state,coeff, modelType)

if (modelType == 1)
    model = xacmes_buildModel_RANKSVM(cur_state,coeff);
end;
