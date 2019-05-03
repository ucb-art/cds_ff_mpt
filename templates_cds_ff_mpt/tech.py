# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING, List, Tuple, Optional, Callable

from bag.layout.tech import TechInfoConfig

from . import config as _config
from .mos.base import MOSTechCDSFFMPT

if TYPE_CHECKING:
    from bag.layout.template import TemplateBase


class TechInfoCDSFFMPT(TechInfoConfig):
    def __init__(self, process_params):
        TechInfoConfig.__init__(self, _config, process_params)

        process_params['layout']['mos_tech_class'] = MOSTechCDSFFMPT(_config, self)
        process_params['layout']['laygo_tech_class'] = None
        process_params['layout']['res_tech_class'] = None

    def get_metal_em_specs(self, layer_name, w, l=-1, vertical=False, **kwargs):
        metal_type = self.get_layer_type(layer_name)
        idc = self._get_metal_idc(metal_type, w, l, vertical, **kwargs)
        irms = self._get_metal_irms(layer_name, w, **kwargs)
        ipeak = float('inf')
        return idc, irms, ipeak

    # noinspection PyUnusedLocal,PyMethodMayBeStatic
    def _get_metal_idc_factor(self, mtype, w, l):
        return 1

    def _get_metal_idc(self, metal_type, w, l, vertical, **kwargs):
        if vertical:
            raise NotImplementedError('Vertical DC current not supported yet')

        inorm, woff = 1.0, 0.0
        idc = inorm * self._get_metal_idc_factor(metal_type, w, l) * (w - woff)
        idc_temp = kwargs.get('dc_temp', self.idc_temp)
        return self.get_idc_scale_factor(idc_temp, metal_type) * idc * 1e-3

    # noinspection PyUnusedLocal
    def _get_metal_irms(self, layer_name, w, **kwargs):
        b = 0.0443
        k, wo, a = 6.0, 0.0, 0.2

        irms_dt = kwargs.get('rms_dt', self.irms_dt)
        irms_ma = (k * irms_dt * (w - wo)**2 * (w - wo + a) / (w - wo + b))**0.5
        return irms_ma * 1e-3

    def get_via_em_specs(self, via_name, bm_layer, tm_layer, via_type='square',
                         bm_dim=(-1, -1), tm_dim=(-1, -1), array=False, **kwargs):
        bm_type = self.get_layer_type(bm_layer)
        tm_type = self.get_layer_type(tm_layer)
        idc = self._get_via_idc(via_name, via_type, bm_type, tm_type, bm_dim,
                                tm_dim, array, **kwargs)
        # via do not have AC current specs
        irms = float('inf')
        ipeak = float('inf')
        return idc, irms, ipeak

    # noinspection PyUnusedLocal
    def _get_via_idc(self, vname, via_type, bm_type, tm_type,
                     bm_dim, tm_dim, array, **kwargs):
        if bm_dim[0] > 0:
            bf = self._get_metal_idc_factor(bm_type, bm_dim[0], bm_dim[1])
        else:
            bf = 1.0

        if tm_dim[0] > 0:
            tf = self._get_metal_idc_factor(tm_type, tm_dim[0], tm_dim[1])
        else:
            tf = 1.0

        factor = min(bf, tf)
        if vname == '1x' or vname == '4':
            if via_type == 'square':
                idc = 0.1 * factor
            elif via_type == 'hrect' or via_type == 'vrect':
                idc = 0.2 * factor
            else:
                # we do not support 2X square via, as it has large
                # spacing rule to square/rectangle vias.
                raise ValueError('Unsupported via type %s' % via_type)
        elif vname == '2x':
            # if via_type == 'square':
                idc = 0.4 * factor
            # else:
            #     raise ValueError('yayUnsupported via type %s' % via_type)
        elif vname == '4x':
                idc = 0.4 * factor
        elif vname == 't':
                idc = 0.4 * factor
        else:
            raise ValueError('Unsupported via name %s and bm_type %s' % (vname, bm_type))

        idc_temp = kwargs.get('dc_temp', self.idc_temp)
        return self.get_idc_scale_factor(idc_temp, bm_type) * idc * 1e-3

    def get_res_em_specs(self, res_type, w, l=-1, **kwargs):
        idc_temp = kwargs.get('dc_temp', self.idc_temp)
        idc_scale = self.get_idc_scale_factor(idc_temp, '', is_res=True)
        idc = 1.0e-3 * w * idc_scale

        irms_dt = kwargs.get('rms_dt', self.irms_dt)
        irms = 1e-3 * (0.02 * irms_dt * w * (w + 0.5)) ** 0.5

        ipeak = 5e-3 * 2 * w
        return idc, irms, ipeak

    def add_cell_boundary(self, template, box):
        pass

    def draw_device_blockage(self, template):
        # type: (TemplateBase) -> None
        pass

    def get_via_arr_enc(self, vname, vtype, mtype, mw_unit, is_bot):
        # type: (...) -> Tuple[Optional[List[Tuple[int, int]]], Optional[Callable[[int, int], bool]]]
        return None, None
