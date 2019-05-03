# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING, Dict, Any, List, Tuple, Union, Optional

from itertools import chain, repeat

from bag.math import lcm
from bag.layout.util import BBox
from bag.layout.template import TemplateBase
from bag.layout.routing import WireArray, TrackID

from abs_templates_ec.analog_mos.finfet import MOSTechFinfetBase

if TYPE_CHECKING:
    from bag.layout.tech import TechInfoConfig


class MOSTechCDSFFMPT(MOSTechFinfetBase):

    def __init__(self, config, tech_info):
        # type: (Dict[str, Any], TechInfoConfig) -> None
        MOSTechFinfetBase.__init__(self, config, tech_info)

    def is_planar_substrate(self, lch_unit, **kwargs):
        # return True to draw planar style substrate
        return False

    def get_conn_yloc_info(self, lch_unit, od_y, md_y, is_sub):
        # type: (int, Tuple[int, int], Tuple[int, int], bool) -> Dict[str, Any]

        mos_constants = self.get_mos_tech_constants(lch_unit)
        mp_h = mos_constants['mp_h']
        mp_h_sub = mos_constants['mp_h_sub']
        mp_md_sp_sub = mos_constants['mp_md_sp_sub']
        ds_m2_sp = mos_constants['ds_m2_sp']

        md_yb, md_yt = md_y
        od_yc = (od_y[0] + od_y[1]) // 2

        # compute gate/drain connection parameters
        g_conn_info = self.get_conn_drc_info(lch_unit, 'g')
        g_m1_h = g_conn_info[1]['min_len']
        g_m1_top_exty = g_conn_info[1]['top_ext']
        g_m1_bot_exty = g_conn_info[1]['bot_ext']
        g_m1_sple = g_conn_info[1]['sp_le']
        g_m2_h = g_conn_info[2]['w']
        g_m3_h = g_conn_info[3]['min_len']
        g_m3_exty = g_conn_info[3]['top_ext']

        d_conn_info = self.get_conn_drc_info(lch_unit, 'd')
        d_m1_bot_exty = d_conn_info[1]['bot_ext']
        d_m2_h = d_conn_info[2]['w']
        d_m3_h = d_conn_info[3]['min_len']
        d_m1_h = max(md_yt - md_yb, d_conn_info[1]['min_len'],
                     2 * d_m1_bot_exty + d_m2_h + ds_m2_sp)

        if is_sub:
            # update MP Y coordinates, compute M1 upper and lower bound
            # bottom MP
            bot_mp_yt = md_yb - mp_md_sp_sub
            bot_mp_yb = bot_mp_yt - mp_h_sub
            bot_mp_yc = (bot_mp_yt + bot_mp_yb) // 2
            g_m1_yb = bot_mp_yc - g_m1_top_exty
            # top MP
            top_mp_yb = md_yt + mp_md_sp_sub
            top_mp_yt = top_mp_yb + mp_h_sub
            top_mp_yc = (top_mp_yb + top_mp_yt) // 2
            g_m1_yt = top_mp_yc + g_m1_top_exty

            m1_y = (g_m1_yb, g_m1_yt)
            m2_y = (od_yc - d_m2_h // 2, od_yc + d_m2_h // 2)
            m3_y = (od_yc - d_m3_h // 2, od_yc + d_m3_h // 2)
            d_y_list = [m1_y, m2_y, m3_y]
            return dict(
                mp_y_list=[(bot_mp_yb, bot_mp_yt), (top_mp_yb, top_mp_yt)],
                g_y_list=[],
                d_y_list=d_y_list,
                s_y_list=d_y_list,
            )
        else:
            d_m1_yb = (md_yb + md_yt - d_m1_h) // 2
            d_m1_yt = d_m1_yb + d_m1_h
            # update gate location
            g_m1_yt = d_m1_yb - g_m1_sple
            g_m1_yb = g_m1_yt - g_m1_h
            g_m2_yc = g_m1_yt - g_m1_bot_exty
            mp_yc = g_m1_yt - g_m1_top_exty
            mp_yt = mp_yc + mp_h // 2
            mp_yb = mp_yt - mp_h

            # compute gate M2/M3 locations
            g_m2_yt = g_m2_yc + g_m2_h // 2
            g_m2_yb = g_m2_yt - g_m2_h
            g_m3_yt = g_m2_yc + g_m3_exty
            g_m3_yb = g_m3_yt - g_m3_h

            # compute source wire location
            s_m2_yc = d_m1_yb + d_m1_bot_exty
            s_m2_yb = s_m2_yc - d_m2_h // 2
            s_m2_yt = s_m2_yb + d_m2_h
            s_m3_yb = s_m2_yc - d_m3_h // 2
            s_m3_yt = s_m3_yb + d_m3_h
            # compute drain wire location
            d_m2_yb = s_m2_yt + ds_m2_sp
            d_m2_yt = d_m2_yb + d_m2_h
            d_m2_yc = (d_m2_yb + d_m2_yt) // 2
            d_m3_yb = d_m2_yc - d_m3_h // 2
            d_m3_yt = d_m3_yb + d_m3_h

            return dict(
                mp_y_list=[(mp_yb, mp_yt)],
                g_y_list=[(g_m1_yb, g_m1_yt), (g_m2_yb, g_m2_yt), (g_m3_yb, g_m3_yt)],
                d_y_list=[(d_m1_yb, d_m1_yt), (d_m2_yb, d_m2_yt), (d_m3_yb, d_m3_yt)],
                s_y_list=[(d_m1_yb, d_m1_yt), (s_m2_yb, s_m2_yt), (s_m3_yb, s_m3_yt)],
            )

    def get_mos_yloc_info(self, lch_unit, w, **kwargs):
        # type: (int, int, **kwargs) -> Dict[str, Any]

        mos_constants = self.get_mos_tech_constants(lch_unit)
        fin_p = mos_constants['mos_pitch']
        od_spy = mos_constants['od_spy']
        mp_cpo_sp = mos_constants['mp_cpo_sp']
        mp_h = mos_constants['mp_h']
        mp_spy = mos_constants['mp_spy']
        cpo_od_sp = mos_constants['cpo_od_sp']
        cpo_h = mos_constants['cpo_h']
        md_h_min = mos_constants['md_h_min']
        md_od_exty = mos_constants['md_od_exty']
        md_spy = mos_constants['md_spy']
        ds_m2_sp = mos_constants['ds_m2_sp']

        od_h = self.get_od_h(lch_unit, w)
        md_h = max(od_h + 2 * md_od_exty, md_h_min)

        # compute gate/drain connection parameters
        g_conn_info = self.get_conn_drc_info(lch_unit, 'g')
        g_m1_top_exty = g_conn_info[1]['top_ext']
        g_m1_sple = g_conn_info[1]['sp_le']

        d_conn_info = self.get_conn_drc_info(lch_unit, 'd')
        d_m3_sple = d_conn_info[3]['sp_le']
        d_m1_bot_exty = d_conn_info[1]['bot_ext']
        d_m2_h = d_conn_info[2]['w']
        d_m1_h = max(md_h, d_conn_info[1]['min_len'],
                     2 * d_m1_bot_exty + d_m2_h + ds_m2_sp)

        # place bottom CPO, compute gate/OD locations
        blk_yb = 0
        cpo_bot_yb = blk_yb - cpo_h // 2
        cpo_bot_yt = cpo_bot_yb + cpo_h
        # get gate via/M1 location
        mp_yb = max(mp_spy // 2, cpo_bot_yt + mp_cpo_sp)
        mp_yt = mp_yb + mp_h
        mp_yc = (mp_yt + mp_yb) // 2
        g_m1_yt = mp_yc + g_m1_top_exty
        # get OD location, round to fin grid.
        od_yb = g_m1_yt + g_m1_sple + d_m1_h // 2 - od_h // 2
        od_yb = self.snap_od_edge(lch_unit, od_yb, False, True)
        od_yt = od_yb + od_h
        od_yc = (od_yb + od_yt) // 2
        md_yb = od_yc - md_h // 2
        md_yt = md_yb + md_h
        # compute top CPO location.
        blk_yt = od_yt + cpo_od_sp + cpo_h // 2
        blk_yt = -(-blk_yt // fin_p) * fin_p

        conn_info = self.get_conn_yloc_info(lch_unit, (od_yb, od_yt), (md_yb, md_yt), False)
        d_m1_yt = conn_info['d_y_list'][0][1]
        d_m3_yb, d_m3_yt = conn_info['d_y_list'][2]
        g_m1_yb = conn_info['g_y_list'][0][0]
        g_m3_yb, g_m3_yt = conn_info['g_y_list'][2]
        return dict(
            blk=(blk_yb, blk_yt),
            po=(blk_yb + cpo_h // 2, blk_yt - cpo_h // 2),
            od=(od_yb, od_yt),
            md=(md_yb, md_yt),
            top_margins=dict(
                od=(blk_yt - od_yt, od_spy),
                md=(blk_yt - md_yt, md_spy),
                m1=(blk_yt - d_m1_yt, g_m1_sple),
                m3=(blk_yt - d_m3_yt, d_m3_sple)
            ),
            bot_margins=dict(
                od=(od_yb - blk_yb, od_spy),
                md=(md_yb - blk_yb, md_spy),
                m1=(g_m1_yb - blk_yb, g_m1_sple),
                m3=(g_m3_yb - blk_yb, d_m3_sple),
            ),
            fill_info={},
            g_conn_y=(g_m3_yb, g_m3_yt),
            d_conn_y=(d_m3_yb, d_m3_yt),
        )

    def get_sub_yloc_info(self, lch_unit, w, **kwargs):
        # type: (int, int, **kwargs) -> Dict[str, Any]
        """Get substrate layout information.

        Layout is quite simple.  We use M0PO to short adjacent S/D together, so dummies can be
        connected using only M2 or below.

        Strategy:

        1. Find bottom M0_PO and bottom OD coordinates from spacing rules.
        #. Find template top coordinate by enforcing symmetry around OD center.
        #. Round up template height to blk_pitch, then recenter OD.
        #. make sure MD/M1 are centered on OD.
        """
        blk_pitch = kwargs['blk_pitch']

        mos_constants = self.get_mos_tech_constants(lch_unit)
        fin_p = mos_constants['mos_pitch']
        od_spy = mos_constants['od_spy']
        mp_cpo_sp_sub = mos_constants['mp_cpo_sp_sub']
        mp_h_sub = mos_constants['mp_h_sub']
        mp_spy_sub = mos_constants['mp_spy_sub']
        mp_md_sp_sub = mos_constants['mp_md_sp_sub']
        cpo_h = mos_constants['cpo_h']
        md_h_min = mos_constants['md_h_min']
        md_od_exty = mos_constants['md_od_exty']
        md_spy = mos_constants['md_spy']

        # compute gate/drain connection parameters
        g_conn_info = self.get_conn_drc_info(lch_unit, 'g')
        g_m1_sple = g_conn_info[1]['sp_le']

        d_conn_info = self.get_conn_drc_info(lch_unit, 'd')
        d_m3_sple = d_conn_info[3]['sp_le']

        od_h = self.get_od_h(lch_unit, w)
        md_h = max(od_h + 2 * md_od_exty, md_h_min)

        # figure out Y coordinate of bottom CPO
        cpo_bot_yt = cpo_h // 2
        # find bottom M0_PO coordinate
        mp_yb = max(mp_spy_sub // 2, cpo_bot_yt + mp_cpo_sp_sub)
        mp_yt = mp_yb + mp_h_sub
        # find first-pass OD coordinate
        od_yb = mp_yt + mp_md_sp_sub + md_h // 2 - od_h // 2
        od_yb = self.snap_od_edge(lch_unit, od_yb, False, True)
        od_yt = od_yb + od_h
        cpo_top_yc = od_yt + od_yb
        # fix substrate height quantization, then recenter OD location
        blk_pitch = lcm([blk_pitch, fin_p])
        blk_h = -(-cpo_top_yc // blk_pitch) * blk_pitch
        cpo_top_yc = blk_h
        od_yb = blk_h // 2 - od_h // 2
        od_yb = self.snap_od_edge(lch_unit, od_yb, False, True)
        od_yt = od_yb + od_h
        od_yc = (od_yb + od_yt) // 2
        # find MD Y coordinates
        md_yb = od_yc - md_h // 2
        md_yt = md_yb + md_h
        blk_yb, blk_yt = 0, cpo_top_yc

        conn_info = self.get_conn_yloc_info(lch_unit, (od_yb, od_yt), (md_yb, md_yt), True)
        m1_yb, m1_yt = conn_info['d_y_list'][0]
        m3_yb, m3_yt = conn_info['d_y_list'][2]
        return dict(
            blk=(blk_yb, blk_yt),
            po=(blk_yb + cpo_h // 2, blk_yt - cpo_h // 2),
            od=(od_yb, od_yt),
            md=(md_yb, md_yt),
            top_margins=dict(
                od=(blk_yt - od_yt, od_spy),
                md=(blk_yt - md_yt, md_spy),
                m1=(blk_yt - m1_yt, g_m1_sple),
                m3=(blk_yt - m3_yt, d_m3_sple)
            ),
            bot_margins=dict(
                od=(od_yb - blk_yb, od_spy),
                md=(md_yb - blk_yb, md_spy),
                m1=(m1_yb - blk_yb, g_m1_sple),
                m3=(m3_yb - blk_yb, d_m3_sple),
            ),
            fill_info={},
            g_conn_y=(m3_yb, m3_yt),
            d_conn_y=(m3_yb, m3_yt),
        )

    def up_one_layer(self, template, cur_lay, cur_y, via_dim, via_sp, via_ble, via_tle,
                     via_x_list, prev_info, conn_drc_info):
        """A helper method that draws vias to connect to upper layer."""

        res = self.res
        lay_name_table = self.config['layer_name']
        via_id_table = self.config['via_id']

        prev_yb, prev_yt, prev_dir, prev_w, prev_lay_name = prev_info
        cur_yb, cur_yt = cur_y
        via_w, via_h = via_dim

        drc_info = conn_drc_info[cur_lay]
        cur_w = drc_info['w']
        cur_dir = drc_info['direction']
        cur_min_len = drc_info['min_len']
        cur_lay_name = lay_name_table[cur_lay]
        via_id = via_id_table[(prev_lay_name, cur_lay_name)]

        # get via Y coord, via enclosures, number of vias, and metal X interval (if horizontal)
        cur_xl = cur_xr = None
        bot_ency, top_ency = via_ble, via_tle
        bot_encx, top_encx = (prev_w - via_w) // 2, (cur_w - via_w) // 2
        if prev_dir == cur_dir:
            # must be both vertical
            arr_yb = max(prev_yb + via_ble, cur_yb + via_tle)
            arr_yt = min(prev_yt - via_ble, cur_yt - via_tle)
            num_rows = (arr_yt - arr_yb + via_sp) // (via_h + via_sp)
            via_yc = (arr_yt + arr_yb) // 2
        else:
            num_rows = 1
            if cur_dir == 'x':
                via_yc = (cur_yb + cur_yt) // 2
                top_encx, top_ency = via_tle, (cur_w - via_h) // 2
                extx = via_w // 2 + top_encx
                # get metal X interval
                cur_xl, cur_xr = via_x_list[0] - extx, via_x_list[-1] + extx
                if cur_min_len > cur_xr - cur_xl:
                    cur_xl = (cur_xr + cur_xl - cur_min_len) // 2
                    cur_xr = cur_xl + cur_min_len
            else:
                via_yc = (prev_yb + prev_yt) // 2
                bot_encx, bot_ency = via_ble, (prev_w - via_h) // 2

        # draw vias and wire(s)
        enc1 = [bot_encx, bot_encx, bot_ency, bot_ency]
        enc2 = [top_encx, top_encx, top_ency, top_ency]
        for via_xc in via_x_list:
            template.add_via_primitive(via_id, [via_xc, via_yc], num_rows=num_rows, sp_rows=via_sp,
                                       enc1=enc1, enc2=enc2, cut_width=via_w, cut_height=via_h,
                                       unit_mode=True)
            if cur_dir == 'y':
                template.add_rect(cur_lay_name, BBox(via_xc - cur_w // 2, cur_yb, via_xc + cur_w // 2, cur_yt,
                                                     res, unit_mode=True))
        if cur_dir == 'x':
            template.add_rect(cur_lay_name, BBox(cur_xl, cur_yb, cur_xr, cur_yt, res, unit_mode=True))

        # setup next iteration
        return cur_yb, cur_yt, cur_dir, cur_w, cur_lay_name

    def draw_ds_connection(self,  # type: MOSTechCDSFFMPT
                           template,  # type: TemplateBase
                           lch_unit,  # type: int
                           fg,  # type: int
                           wire_pitch,  # type: int
                           xc,  # type: int
                           od_y,  # type: Tuple[int, int]
                           md_y,  # type: Tuple[int, int]
                           dum_x_list,  # type: List[int]
                           conn_x_list,  # type: List[int]
                           align_gate,  # type: bool
                           wire_dir,  # type: int
                           ds_code,  # type: int
                           **kwargs
                           ):
        # type: (...) -> Tuple[List[WireArray], List[WireArray]]

        is_dum = kwargs.get('is_dum', False)

        res = self.res
        mos_lay_table = self.config['mos_layer_table']

        mos_constants = self.get_mos_tech_constants(lch_unit)
        md_w = mos_constants['md_w']
        bot_layer = mos_constants['d_bot_layer']
        via_info = mos_constants['d_via']

        is_sub = (ds_code == 3)
        conn_yloc_info = self.get_conn_yloc_info(lch_unit, od_y, md_y, is_sub)
        conn_drc_info = self.get_conn_drc_info(lch_unit, 'd')

        dum_layer = self.get_dum_conn_layer()
        mos_layer = self.get_mos_conn_layer()
        dum_warrs, conn_warrs = [], []

        # figure out via X coordinates
        if is_sub:
            via_x_list = list(range(xc, xc + (fg + 1) * wire_pitch, wire_pitch))
            conn_y_list = conn_yloc_info['d_y_list']
        else:
            if ds_code == 1:
                via_x_list = list(range(xc, xc + (fg + 1) * wire_pitch, 2 * wire_pitch))
            else:
                via_x_list = list(range(xc + wire_pitch, xc + (fg + 1) * wire_pitch, 2 * wire_pitch))
            if align_gate:
                conn_y_list = conn_yloc_info['d_y_list']
            else:
                conn_y_list = conn_yloc_info['s_y_list']

        # connect from OD up to M3
        stop_layer = dum_layer if is_dum else mos_layer
        lay_list = range(bot_layer, stop_layer + 1)
        prev_info = md_y[0], md_y[1], 'y', md_w, mos_lay_table['MD']
        for cur_lay, cur_y, via_dim, via_sp, via_ble, via_tle in \
                zip(lay_list, conn_y_list, via_info['dim'], via_info['sp'],
                    via_info['bot_enc_le'], via_info['top_enc_le']):
            prev_info = self.up_one_layer(template, cur_lay, cur_y, via_dim, via_sp, via_ble, via_tle,
                                          via_x_list, prev_info, conn_drc_info)

        # add WireArrays
        if stop_layer >= dum_layer:
            cur_yb, cur_yt = conn_y_list[dum_layer - bot_layer]
            for conn_xc in dum_x_list:
                tidx = template.grid.coord_to_track(dum_layer, conn_xc, unit_mode=True)
                dum_warrs.append(WireArray(TrackID(dum_layer, tidx), cur_yb * res, cur_yt * res, res))
        if stop_layer >= mos_layer:
            cur_yb, cur_yt = conn_y_list[mos_layer - bot_layer]
            for conn_xc in conn_x_list:
                tidx = template.grid.coord_to_track(mos_layer, conn_xc, unit_mode=True)
                conn_warrs.append(WireArray(TrackID(mos_layer, tidx), cur_yb * res, cur_yt * res, res))

        return dum_warrs, conn_warrs

    def draw_g_connection(self,  # type: MOSTechCDSFFMPT
                          template,  # type: TemplateBase
                          lch_unit,  # type: int
                          fg,  # type: int
                          sd_pitch,  # type: int
                          xc,  # type: int
                          od_y,  # type: Tuple[int, int]
                          md_y,  # type: Tuple[int, int]
                          conn_x_list,  # type: List[int]
                          is_sub=False,  # type: bool
                          **kwargs
                          ):
        # type: (...) -> List[WireArray]

        is_dum = kwargs.get('is_dum', False)
        res = self.res
        mos_lay_table = self.config['mos_layer_table']
        lay_name_table = self.config['layer_name']
        via_id_table = self.config['via_id']

        mos_constants = self.get_mos_tech_constants(lch_unit)
        mp_po_ovl_constants = mos_constants['mp_po_ovl_constants']
        mp_po_ovl_constants_sub = mos_constants['mp_po_ovl_constants_sub']
        mp_h = mos_constants['mp_h']
        mp_h_sub = mos_constants['mp_h_sub']
        via_info = mos_constants['g_via']

        conn_yloc_info = self.get_conn_yloc_info(lch_unit, od_y, md_y, is_sub)
        conn_drc_info = self.get_conn_drc_info(lch_unit, 'g')

        conn_warrs = []

        mp_lay = mos_lay_table['MP']
        m1_w = conn_drc_info[1]['w']
        mp_y_list = conn_yloc_info['mp_y_list']
        v0_id = via_id_table[(mos_lay_table['MP'], lay_name_table[1])]
        if is_sub:
            mp_po_ovl = mp_po_ovl_constants_sub[0] + lch_unit * mp_po_ovl_constants_sub[1]
            # connect gate to M1 only
            m1_yb, m1_yt = conn_yloc_info['d_y_list'][0]
            via_w, via_h = via_info['dim'][0]
            bot_encx = via_info['bot_enc_le'][0]
            top_encx = (m1_w - via_w) // 2
            bot_ency = (mp_h_sub - via_h) // 2
            top_ency = via_info['top_enc_le'][0]

            mp_dx = sd_pitch // 2 - lch_unit // 2 + mp_po_ovl
            enc1 = [bot_encx, bot_encx, bot_ency, bot_ency]
            enc2 = [top_encx, top_encx, top_ency, top_ency]
            for via_xc in range(xc, xc + (fg + 1) * sd_pitch, 2 * sd_pitch):
                mp_xl = via_xc - mp_dx
                mp_xr = via_xc + mp_dx
                for mp_yb, mp_yt in mp_y_list:
                    template.add_rect(mp_lay, BBox(mp_xl, mp_yb, mp_xr, mp_yt, res, unit_mode=True))
                    mp_yc = (mp_yb + mp_yt) // 2
                    template.add_via_primitive(v0_id, [via_xc, mp_yc], enc1=enc1, enc2=enc2,
                                               cut_width=via_w, cut_height=via_h, unit_mode=True)
                template.add_rect('M1', BBox(via_xc - m1_w // 2, m1_yb, via_xc + m1_w // 2, m1_yt, res,
                                             unit_mode=True))
        else:
            mp_po_ovl = mp_po_ovl_constants[0] + lch_unit * mp_po_ovl_constants[1]

            if fg % 2 == 0:
                gate_fg_list = [2] * (fg // 2)
            else:
                if fg == 1:
                    raise ValueError('cannot connect 1 finger transistor')
                if fg <= 5:
                    gate_fg_list = [fg]
                elif (fg - 3) % 4 == 0:
                    num_mp_half = (fg - 3) // 2
                    gate_fg_list = list(chain(repeat(2, num_mp_half // 2), [3], repeat(2, num_mp_half // 2)))
                else:
                    num_mp_half = (fg - 5) // 2
                    gate_fg_list = list(chain(repeat(2, num_mp_half // 2), [5], repeat(2, num_mp_half // 2)))

            # connect gate to M1.
            tot_fg = 0
            mp_yb, mp_yt = mp_y_list[0]
            m1_yb, m1_yt = conn_yloc_info['g_y_list'][0]
            via_w, via_h = via_info['dim'][0]
            bot_encx = via_info['bot_enc_le'][0]
            top_encx = (m1_w - via_w) // 2
            bot_ency = (mp_h - via_h) // 2
            top_ency = via_info['top_enc_le'][0]
            enc1 = [bot_encx, bot_encx, bot_ency, bot_ency]
            enc2 = [top_encx, top_encx, top_ency, top_ency]
            via_yc = (mp_yb + mp_yt) // 2
            via_x_list = []
            for num_fg in gate_fg_list:
                via_xoff = xc + (tot_fg + 1) * sd_pitch
                cur_xc = xc + tot_fg * sd_pitch + num_fg * sd_pitch // 2
                # draw MP
                mp_w = (num_fg - 1) * sd_pitch - lch_unit + 2 * mp_po_ovl
                mp_xl = cur_xc - mp_w // 2
                mp_xr = mp_xl + mp_w
                template.add_rect(mp_lay, BBox(mp_xl, mp_yb, mp_xr, mp_yt, res, unit_mode=True))
                # draw V0, M1
                for via_xc in range(via_xoff, via_xoff + (num_fg - 1) * sd_pitch, sd_pitch):
                    cur_tidx = template.grid.coord_to_track(1, via_xc, unit_mode=True)
                    template.add_via_primitive(v0_id, [via_xc, via_yc], enc1=enc1, enc2=enc2,
                                               cut_width=via_w, cut_height=via_h, unit_mode=True)
                    template.add_wires(1, cur_tidx, m1_yb, m1_yt, unit_mode=True)
                    via_x_list.append(via_xc)
                tot_fg += num_fg

            # connect from M1 up to M3 if not dummy gate connection
            if not is_dum:
                conn_y_list = conn_yloc_info['g_y_list'][1:]
                lay_list = range(2, 2 + len(conn_y_list))
                prev_info = m1_yb, m1_yt, 'y', m1_w, lay_name_table[1]
                for cur_lay, cur_y, via_dim, via_sp, via_ble, via_tle in \
                        zip(lay_list, conn_y_list, via_info['dim'][1:], via_info['sp'][1:],
                            via_info['bot_enc_le'][1:], via_info['top_enc_le'][1:]):
                    prev_info = self.up_one_layer(template, cur_lay, cur_y, via_dim, via_sp, via_ble, via_tle,
                                                  via_x_list, prev_info, conn_drc_info)
                    via_x_list = conn_x_list

                # add ports
                mos_layer = self.get_mos_conn_layer()
                cur_yb, cur_yt = conn_y_list[-1]
                for conn_xc in conn_x_list:
                    tidx = template.grid.coord_to_track(mos_layer, conn_xc, unit_mode=True)
                    conn_warrs.append(WireArray(TrackID(mos_layer, tidx), cur_yb * res, cur_yt * res, res))

        return conn_warrs

    def draw_substrate_connection(self,  # type: MOSTechFinfetBase
                                  template,  # type: TemplateBase
                                  layout_info,  # type: Dict[str, Any]
                                  port_tracks,  # type: List[Union[float, int]]
                                  dum_tracks,  # type: List[Union[float, int]]
                                  exc_tracks,  # type: List[Union[float, int]]
                                  dummy_only,  # type: bool
                                  is_laygo,  # type: bool
                                  is_guardring,  # type: bool
                                  options,  # type: Dict[str, Any]
                                  ):
        # type: (...) -> bool

        sub_parity = options.get('sub_parity', 0)

        lch_unit = layout_info['lch_unit']
        row_info_list = layout_info['row_info_list']

        mos_constants = self.get_mos_tech_constants(lch_unit)
        sd_pitch = mos_constants['sd_pitch']

        sd_pitch2 = sd_pitch // 2

        exc_set = set(int(2 * v + 1) for v in exc_tracks)

        has_od = False
        for row_info in row_info_list:
            od_y = row_info.od_y
            md_y = row_info.md_y
            if od_y[1] > od_y[0]:
                has_od = True

                # find current port name
                od_start, od_stop = row_info.od_x_list[0]
                fg = od_stop - od_start
                xshift = od_start * sd_pitch
                sub_type = row_info.od_type[1]
                port_name = 'VDD' if sub_type == 'ntap' else 'VSS'

                # find X locations of M1/M3.
                if dummy_only:
                    # find X locations to draw vias
                    dum_x_list = [sd_pitch2 * int(2 * v + 1) for v in dum_tracks]
                    conn_x_list = []
                else:
                    # first, figure out port/dummy tracks
                    # To lower parasitics, we try to draw only as many dummy tracks as necessary.
                    # Also, every port track is also a dummy track (because some technology
                    # there's no horizontal short).  With these constraints, our track selection
                    # algorithm is as follows:
                    # 1. for every dummy track, if its not adjacent to any port tracks, add it to
                    #    port tracks (this improves dummy connection resistance to supply).
                    # 2. Try to add as many unused tracks to port tracks as possible, while making
                    #    sure we don't end up with adjacent port tracks.  This improves substrate
                    #    connection resistance to supply.

                    # use half track indices so we won't have rounding errors.
                    dum_htr_set = set((int(2 * v + 1) for v in dum_tracks))
                    conn_htr_set = set((int(2 * v + 1) for v in port_tracks))
                    # add as many dummy tracks as possible to port tracks
                    for d in dum_htr_set:
                        if d + 2 not in conn_htr_set and d - 2 not in conn_htr_set:
                            conn_htr_set.add(d)
                    # add as many unused tracks as possible to port tracks
                    for htr in range(od_start * 2, 2 * od_stop + 1, 2):
                        if htr + 2 not in conn_htr_set and htr - 2 not in conn_htr_set:
                            conn_htr_set.add(htr)
                    # add all port sets to dummy set
                    dum_htr_set.update(conn_htr_set)
                    # find X coordinates
                    dum_x_list = [sd_pitch2 * v for v in sorted(dum_htr_set)]
                    conn_x_list = [sd_pitch2 * v for v in sorted(conn_htr_set)]

                ds_code = 4 if is_guardring else 3
                dum_warrs, port_warrs = self.draw_ds_connection(template, lch_unit, fg, sd_pitch,
                                                                xshift, od_y, md_y,
                                                                dum_x_list, conn_x_list, True, 1,
                                                                ds_code,
                                                                ud_parity=sub_parity)
                template.add_pin(port_name, dum_warrs, show=False)
                template.add_pin(port_name, port_warrs, show=False)

                if not is_guardring:
                    self.draw_g_connection(template, lch_unit, fg, sd_pitch, xshift, od_y, md_y,
                                           conn_x_list, is_sub=True)
        return has_od

    def draw_dum_connection_helper(self,
                                   template,  # type: TemplateBase
                                   lch_unit,  # type: int
                                   fg,  # type: int
                                   sd_pitch,  # type: int
                                   xc,  # type: int
                                   od_y,  # type: Tuple[int, int]
                                   md_y,  # type: Tuple[int, int]
                                   ds_x_list,  # type: List[int]
                                   gate_tracks,  # type: List[Union[float, int]]
                                   left_edge,  # type: bool
                                   right_edge,  # type: bool
                                   options,  # type: Dict[str, Any]
                                   ):
        # type: (...) -> List[WireArray]

        res = self.res
        lay_name_table = self.config['layer_name']

        mos_constants = self.get_mos_tech_constants(lch_unit)
        g_m1_dum_h = mos_constants['g_m1_dum_h']

        conn_yloc_info = self.get_conn_yloc_info(lch_unit, od_y, md_y, False)

        m1_lay = lay_name_table[1]
        m1_yb = conn_yloc_info['g_y_list'][0][0]
        m1_yt = conn_yloc_info['d_y_list'][0][1]

        # draw gate/drain/source connection to M1
        self.draw_g_connection(template, lch_unit, fg, sd_pitch, xc, od_y, md_y, [], is_sub=False, is_dum=True)
        self.draw_ds_connection(template, lch_unit, fg, sd_pitch, xc, od_y, md_y, [], [], False, 1, 1, is_dum=True)
        self.draw_ds_connection(template, lch_unit, fg, sd_pitch, xc, od_y, md_y, [], [], True, 1, 2, is_dum=True)

        # short M1 together
        dum_layer = self.get_dum_conn_layer()
        for wire_xc in ds_x_list:
            tidx = template.grid.coord_to_track(dum_layer, wire_xc, unit_mode=True)
            template.add_wires(dum_layer, tidx, m1_yb, m1_yt, unit_mode=True)
        xl = ds_x_list[0]
        xr = ds_x_list[-1]
        if xr > xl:
            template.add_rect(m1_lay, BBox(xl, m1_yb, xr, m1_yb + g_m1_dum_h, res, unit_mode=True))

        # return gate ports
        return [WireArray(TrackID(dum_layer, tidx), m1_yb * res, m1_yt * res, res) for tidx in gate_tracks]

    def draw_decap_connection_helper(self,
                                     template,  # type: TemplateBase
                                     lch_unit,  # type: int
                                     fg,  # type: int
                                     sd_pitch,  # type: int
                                     xc,  # type: int
                                     od_y,  # type: Tuple[int, int]
                                     md_y,  # type: Tuple[int, int]
                                     gate_ext_mode,  # type: int
                                     export_gate,  # type: bool
                                     ):
        # type: (...) -> Tuple[Optional[WireArray], List[WireArray]]
        raise NotImplementedError('Not implemented.')
