function [best_val1, best_val2, best_val3, best_meas] = checkbesthyperFromcoarse_fine()
meas_name = 'snr_dynamic'; %'error_trend_sign_all';
out_dir = 'Root1658g145TF_Quant_df2_1e2_output_coarse_fine/';
files = dir([out_dir '/*.mat']);
par1 = 'gamma';
par2 = 'lambda_w';
par3= 'tau';
% Inspect all results
best_meas = -inf;
best_val1 = -inf;
best_val2 = -inf;
best_val3 = -inf;
for i = 1:length(files)
  % Get the MODEL and METER_INFER_TEST variables from file <i>
  d = load(sprintf('%s/%s', out_dir, files(i).name));
  meas = eval(['d.METER_INFER_TEST{end-1}.last_' meas_name ';']);
  if (meas > best_meas)
    best_meas = meas;

    % Evaluate first parameter
    best_val1 = eval(['d.params.' par1 ';']);

    % Evaluate second parameter?
    if ~isempty(par2)
      best_val2 = eval(['d.params.' par2 ';']);
    end

    % Evaluate third parameter?
    if ~isempty(par3)
      best_val3 = eval(['d.params.' par3 ';']);
    end
  end
end

end
