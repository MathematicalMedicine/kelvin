/* Copyright (C) 2012, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
#include <stdio.h>
#include <stdlib.h>
#include "dcuhre.h"
#include "utils/tpl.h"

int main(int argc, char *argv[]) {
  tpl_node *tn;
  char *regionTPLFormat = "iiffii";
  FILE *file;
  sub_region s;

  if (argc != 2) {
    fprintf (stderr, "Usage: %s <region file name>\n", argv[0]);
    exit (1);
  }
  tn = tpl_map (regionTPLFormat, &s.parent_id, &s.region_level, &s.local_result, &s.local_error, &s.dir, &s.cur_scale);
  if ((file = fopen (argv[1], "r")) == NULL) {
    fprintf (stderr, "File %s cannot be read\n", argv[1]);
    exit (2);
  }
  fclose (file);
  tpl_load (tn, TPL_FILE, argv[1]);
  tpl_unpack (tn, 0);
  tpl_free (tn);
  fprintf (stdout, "parent_id: %d, region_level: %d, local_result: %g, local_error: %g, dir: %d, cur_scale: %d\n",
	   s.parent_id, s.region_level, s.local_result, s.local_error, s.dir, s.cur_scale);
  return 0;
}


