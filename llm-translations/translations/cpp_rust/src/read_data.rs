#![allow(dead_code)]

// Direct port of readData.cpp (TSPLIB-style instance reader). Not called from
// main, exactly as in the original C++ (main.cpp reads the pre-built
// distance matrix via data::load_data instead).

fn next_token(tokens: &[&str], idx: &mut usize) -> String {
    let t = tokens.get(*idx).copied().unwrap_or("");
    *idx += 1;
    t.to_string()
}

pub fn read_data(args: &[String]) -> (i32, Vec<Vec<f64>>) {
    if args.len() < 2 {
        println!("\nFaltando parametros");
        println!(" ./exec [Instancia] ");
        std::process::exit(1);
    }

    if args.len() > 2 {
        println!("\nMuitos parametros");
        println!(" ./exec [Instancia] ");
        std::process::exit(1);
    }

    let instancia = &args[1];

    let content = match std::fs::read_to_string(instancia) {
        Ok(c) => c,
        Err(_) => {
            println!("arquivo nao pode ser aberto");
            std::process::exit(1);
        }
    };

    let tokens: Vec<&str> = content.split_whitespace().collect();
    let mut idx: usize = 0;

    let mut arquivo = String::new();

    while arquivo != "DIMENSION:" && arquivo != "DIMENSION" {
        arquivo = next_token(&tokens, &mut idx);
    }
    if arquivo == "DIMENSION" {
        arquivo = next_token(&tokens, &mut idx);
    }

    let n: i32 = next_token(&tokens, &mut idx).parse().unwrap();

    while arquivo != "EDGE_WEIGHT_TYPE:" && arquivo != "EDGE_WEIGHT_TYPE" {
        arquivo = next_token(&tokens, &mut idx);
    }
    if arquivo == "EDGE_WEIGHT_TYPE" {
        arquivo = next_token(&tokens, &mut idx);
    }

    let ewt = next_token(&tokens, &mut idx);

    let nn = (n + 1) as usize;
    let mut x = vec![0.0f64; nn];
    let mut y = vec![0.0f64; nn];
    let mut dist = vec![vec![0.0f64; nn]; nn];

    if ewt == "EXPLICIT" {
        while arquivo != "EDGE_WEIGHT_FORMAT:" && arquivo != "EDGE_WEIGHT_FORMAT" {
            arquivo = next_token(&tokens, &mut idx);
        }
        if arquivo == "EDGE_WEIGHT_FORMAT" {
            arquivo = next_token(&tokens, &mut idx);
        }

        let ewf = next_token(&tokens, &mut idx);

        if ewf == "FUNCTION" {
            println!("FUNCTION - Nao suportado!");
        } else if ewf == "FULL_MATRIX" {
            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = next_token(&tokens, &mut idx);
            }

            for i in 1..(n as usize + 1) {
                for j in 1..(n as usize + 1) {
                    dist[i][j] = next_token(&tokens, &mut idx).parse().unwrap();
                }
            }
        } else if ewf == "UPPER_ROW" {
            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = next_token(&tokens, &mut idx);
            }

            let mut i = 1usize;
            while i < n as usize {
                let mut j = i + 1;
                while j < n as usize + 1 {
                    let v: f64 = next_token(&tokens, &mut idx).parse().unwrap();
                    dist[i][j] = v;
                    dist[j][i] = v;
                    j += 1;
                }
                i += 1;
            }

            for i in 1..(n as usize + 1) {
                dist[i][i] = 0.0;
            }
        } else if ewf == "LOWER_ROW" {
            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = next_token(&tokens, &mut idx);
            }

            let mut i = 2usize;
            while i < n as usize + 1 {
                let mut j = 1usize;
                while j < i {
                    let v: f64 = next_token(&tokens, &mut idx).parse().unwrap();
                    dist[i][j] = v;
                    dist[j][i] = v;
                    j += 1;
                }
                i += 1;
            }

            for i in 1..(n as usize + 1) {
                dist[i][i] = 0.0;
            }
        } else if ewf == "UPPER_DIAG_ROW" {
            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = next_token(&tokens, &mut idx);
            }

            let mut i = 1usize;
            while i < n as usize + 1 {
                let mut j = i;
                while j < n as usize + 1 {
                    let v: f64 = next_token(&tokens, &mut idx).parse().unwrap();
                    dist[i][j] = v;
                    dist[j][i] = v;
                    j += 1;
                }
                i += 1;
            }
        } else if ewf == "LOWER_DIAG_ROW" {
            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = next_token(&tokens, &mut idx);
            }

            let mut i = 1usize;
            while i < n as usize + 1 {
                let mut j = 1usize;
                while j <= i {
                    let v: f64 = next_token(&tokens, &mut idx).parse().unwrap();
                    dist[i][j] = v;
                    dist[j][i] = v;
                    j += 1;
                }
                i += 1;
            }
        } else if ewf == "UPPER_COL" {
            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = next_token(&tokens, &mut idx);
            }

            let mut j = 2usize;
            while j < n as usize + 1 {
                let mut i = 1usize;
                while i < j {
                    let v: f64 = next_token(&tokens, &mut idx).parse().unwrap();
                    dist[i][j] = v;
                    dist[j][i] = v;
                    i += 1;
                }
                j += 1;
            }

            for i in 1..(n as usize + 1) {
                dist[i][i] = 0.0;
            }
        } else if ewf == "LOWER_COL" {
            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = next_token(&tokens, &mut idx);
            }

            // NOTE: mirrors the original C++ loop verbatim, including its bug
            // (the inner loop increments `j` instead of `i`, matching
            // `for (int i = j+1; i < N+1; j++)` in readData.cpp).
            let mut j = 1usize;
            while j < n as usize {
                let i = j + 1;
                while i < n as usize + 1 {
                    let v: f64 = next_token(&tokens, &mut idx).parse().unwrap();
                    dist[i][j] = v;
                    dist[j][i] = v;
                    j += 1;
                }
                let _ = i;
            }

            for i in 1..(n as usize + 1) {
                dist[i][i] = 0.0;
            }
        } else if ewf == "UPPER_DIAG_COL" {
            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = next_token(&tokens, &mut idx);
            }

            let mut j = 1usize;
            while j < n as usize + 1 {
                let mut i = 1usize;
                while i <= j {
                    let v: f64 = next_token(&tokens, &mut idx).parse().unwrap();
                    dist[i][j] = v;
                    dist[j][i] = v;
                    i += 1;
                }
                j += 1;
            }
        } else if ewf == "LOWER_DIAG_COL" {
            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = next_token(&tokens, &mut idx);
            }

            // NOTE: mirrors the original C++ loop verbatim, including its bug
            // (the inner loop increments `j` instead of `i`, matching
            // `for (int i = j; i < N+1; j++)` in readData.cpp).
            let mut j = 1usize;
            while j < n as usize + 1 {
                let i = j;
                while i < n as usize + 1 {
                    let v: f64 = next_token(&tokens, &mut idx).parse().unwrap();
                    dist[i][j] = v;
                    dist[j][i] = v;
                    j += 1;
                }
                let _ = i;
            }
        }
    } else if ewt == "EUC_2D" {
        while arquivo != "NODE_COORD_SECTION" {
            arquivo = next_token(&tokens, &mut idx);
        }

        for i in 1..(n as usize + 1) {
            let _temp_city = next_token(&tokens, &mut idx);
            x[i] = next_token(&tokens, &mut idx).parse().unwrap();
            y[i] = next_token(&tokens, &mut idx).parse().unwrap();
        }

        for i in 1..(n as usize + 1) {
            for j in 1..(n as usize + 1) {
                dist[i][j] = (calc_dist_euc(&x, &y, i, j) + 0.5).floor();
            }
        }
    } else if ewt == "EUD_3D" {
        println!("EUC_3D - Nao suportado!");
    } else if ewt == "MAX_2D" {
        println!("MAX_2D - Nao suportado!");
    } else if ewt == "MAX_3D" {
        println!("MAX_3D - Nao Suportado!");
    } else if ewt == "MAN_2D" {
        println!("MAN_2D - Nao suportado!");
    } else if ewt == "MAN_3D" {
        println!("MAN_3D - Nao Suportado!");
    } else if ewt == "CEIL_2D" {
        while arquivo != "NODE_COORD_SECTION" {
            arquivo = next_token(&tokens, &mut idx);
        }

        for i in 1..(n as usize + 1) {
            let _temp_city = next_token(&tokens, &mut idx);
            x[i] = next_token(&tokens, &mut idx).parse().unwrap();
            y[i] = next_token(&tokens, &mut idx).parse().unwrap();
        }

        for i in 1..(n as usize + 1) {
            for j in 1..(n as usize + 1) {
                dist[i][j] = calc_dist_euc(&x, &y, i, j).ceil();
            }
        }
    } else if ewt == "GEO" {
        while arquivo != "NODE_COORD_SECTION" {
            arquivo = next_token(&tokens, &mut idx);
        }

        for i in 1..(n as usize + 1) {
            let _temp_city = next_token(&tokens, &mut idx);
            x[i] = next_token(&tokens, &mut idx).parse().unwrap();
            y[i] = next_token(&tokens, &mut idx).parse().unwrap();
        }

        let mut latitude = vec![0.0f64; nn];
        let mut longitude = vec![0.0f64; nn];

        calc_lat_long(&x, &y, n, &mut latitude, &mut longitude);

        for i in 1..(n as usize + 1) {
            for j in 1..(n as usize + 1) {
                dist[i][j] = calc_dist_geo(&latitude, &longitude, i, j);
            }
        }
    } else if ewt == "ATT" {
        while arquivo != "NODE_COORD_SECTION" {
            arquivo = next_token(&tokens, &mut idx);
        }

        let mut temp_x = vec![0i32; nn];
        let mut temp_y = vec![0i32; nn];

        for i in 1..(n as usize + 1) {
            let _temp_city = next_token(&tokens, &mut idx);
            temp_x[i] = next_token(&tokens, &mut idx).parse().unwrap();
            temp_y[i] = next_token(&tokens, &mut idx).parse().unwrap();
            x[i] = temp_x[i] as f64;
            y[i] = temp_y[i] as f64;
        }

        for i in 1..(n as usize + 1) {
            for j in 1..(n as usize + 1) {
                dist[i][j] = calc_dist_att(&x, &y, i, j);
            }
        }
    } else if ewt == "XRAY1" {
        println!("XRAY1 - Nao suportado!");
    } else if ewt == "XRAY2" {
        println!("XRAY2 - Nao suportado!");
    } else if ewt == "SPECIAL" {
        println!("SPECIAL - Nao suportado!");
    }

    (n, dist)
}

fn calc_dist_euc(x: &[f64], y: &[f64], i: usize, j: usize) -> f64 {
    ((x[i] - x[j]).powi(2) + (y[i] - y[j]).powi(2)).sqrt()
}

fn calc_dist_att(x: &[f64], y: &[f64], i: usize, j: usize) -> f64 {
    let rij = (((x[i] - x[j]).powi(2) + (y[i] - y[j]).powi(2)) / 10.0).sqrt();
    let tij = (rij + 0.5).floor();

    if tij < rij {
        tij + 1.0
    } else {
        tij
    }
}

fn calc_lat_long(x: &[f64], y: &[f64], n: i32, latit: &mut [f64], longit: &mut [f64]) {
    let pi = 3.141592_f64;

    for i in 1..(n as usize + 1) {
        let deg = x[i].trunc();
        let min = x[i] - deg;
        latit[i] = pi * (deg + 5.0 * min / 3.0) / 180.0;
    }

    for i in 1..(n as usize + 1) {
        let deg = y[i].trunc();
        let min = y[i] - deg;
        longit[i] = pi * (deg + 5.0 * min / 3.0) / 180.0;
    }
}

fn calc_dist_geo(latit: &[f64], longit: &[f64], i: usize, j: usize) -> f64 {
    let rrr = 6378.388_f64;

    let q1 = (longit[i] - longit[j]).cos();
    let q2 = (latit[i] - latit[j]).cos();
    let q3 = (latit[i] + latit[j]).cos();

    ((rrr * (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)).acos()) + 1.0).trunc()
}
