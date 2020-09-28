use crate::types::frame::FrameElement;
use crate::{ANALYSIS_DATA, ELEMENT_DATA};
use na::{MatrixN, U12};
use std::collections::HashMap;

/**
 * Returns the world displacement vectors for all free degrees of freedom.
 *
 * Solves equations of the form `F = k Δ` where
 * * `F` is the world force vector
 * * `k` is the world stiffness matrix
 * * `Δ` is the world displacement vector
 *
 * # Arguments
 *
 * `forces` - forces acting at the free degrees of freedom
 * `stiffness` - world stiffness matrix for the free degrees of freedom
 */
pub fn solve_for_deflections(
    stiffness: &mut sparse21::Matrix,
    forces: Vec<f64>,
) -> Result<Vec<f64>, String> {
    let solution: Result<Vec<f64>, &str> = stiffness.solve(forces);

    match solution {
        Ok(result) => Ok(result),
        Err(e) => Err(e.to_string()),
    }
}

/**
 * Returns the non-zero entries comprising the assembled stiffness matrix in world coordinates.
 *
 * The individual stiffness contributions of each finite element stiffness matrix are summed
 * and combined in correspondence to the world degrees of freedom.
 */
pub fn assemble_world_stiffness_matrix() -> Result<HashMap<(usize, usize), f64>, String> {
    let mut assembled = HashMap::<(usize, usize), f64>::new();

    for frame in (&ELEMENT_DATA.frames).values() {
        merge_stiffness_matrix_at_frame_dofs(&mut assembled, &frame)?;
    }

    Ok(assembled)
}

fn merge_stiffness_matrix_at_frame_dofs(
    assembled: &mut HashMap<(usize, usize), f64>,
    frame: &FrameElement,
) -> Result<(), String> {
    let stiffness: MatrixN<f64, U12> = match ANALYSIS_DATA.frame_stiffnesses.get(&frame.id) {
        Some(x) => x.world,
        None => return Err(format!("Failed to locate frame with id = {}", &frame.id)),
    };

    let start_dofs = match frame.start_node() {
        Some(x) => x.degrees_of_freedom,
        None => return Err(format!("Failed to find node id = {}", &frame.start_node_id)),
    };
    let end_dofs = match frame.end_node() {
        Some(x) => x.degrees_of_freedom,
        None => return Err(format!("Failed to find node id = {}", &frame.end_node_id)),
    };

    for i in 1..6 {
        for j in 1..12 {
            if j < 6 {
                merge_stiffness_matrix_at_dof(
                    assembled,
                    start_dofs[i],
                    start_dofs[j],
                    stiffness[(i, j)],
                )
            } else {
                merge_stiffness_matrix_at_dof(
                    assembled,
                    start_dofs[i],
                    end_dofs[j],
                    stiffness[(i, j)],
                )
            }
        }
    }

    for i in 6..12 {
        for j in 1..12 {
            if j < 6 {
                merge_stiffness_matrix_at_dof(
                    assembled,
                    end_dofs[i],
                    start_dofs[j],
                    stiffness[(i, j)],
                )
            } else {
                merge_stiffness_matrix_at_dof(
                    assembled,
                    end_dofs[i],
                    end_dofs[j],
                    stiffness[(i, j)],
                )
            }
        }
    }

    Ok(())
}

fn merge_stiffness_matrix_at_dof(
    assembled: &mut HashMap<(usize, usize), f64>,
    row: usize,
    column: usize,
    value: f64,
) {
    if let Some(val) = assembled.get_mut(&(row, column)) {
        *val += value;
    } else {
        assembled.insert((row, column), value);
    }
}
