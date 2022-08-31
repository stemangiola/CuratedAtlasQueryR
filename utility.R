
clean_cell_type = function(x){
	x |>
		str_replace_all(" ", "_") |>
		str_replace_all("\\.", "_") |>
		str_replace_all(",", "_") |>
		str_replace_all("\\(", "_") |>
		str_replace_all("\\)", "_") |>
		str_replace_all("/", "_")
}