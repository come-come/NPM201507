0802mor: 
	df = data[(data[data.columns[window]] >= (t - 0.05 - 0.00001)) & (data[data.columns[window]] <= (t + 0.00001))]
	有误差
0802aft:
	df = data[(data[data.columns[window]] >= (t - 0.05 + 0.00001)) & (data[data.columns[window]] <= (t + 0.00001))]