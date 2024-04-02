/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_COLORS_H
#define MECHSYS_COLORS_H

// Std Lib
#include <map>
#include <cstdlib> // for rand

// MechSys
#include <mechsys/util/string.h>
#include <mechsys/linalg/matvec.h>

namespace Colors
{

typedef std::map<String, Vec3_t> CLR_t;

CLR_t CLR;

Vec3_t const & Get (String const & Name)
{
    Colors::CLR_t::const_iterator p = Colors::CLR.find(Name);
    if (p==Colors::CLR.end()) throw new Fatal("Colors::Get: Could not find color named %s in Colors::CLR map",Name.CStr());
    return p->second;
}

Vec3_t const & Get (char const * Name)
{
    Colors::CLR_t::const_iterator p = Colors::CLR.find(Name);
    if (p==Colors::CLR.end()) throw new Fatal("Colors::Get: Could not find color named %s in Colors::CLR map",Name);
    return p->second;
}

void GetRandom (size_t Num, Array<String> & Clrs)
{
    size_t num   = (Num<CLR.size() ? Num : CLR.size());
    size_t maxit = 10*CLR.size();
    size_t iter  = 0;
    while (iter<maxit)
    {
        int count = 0;
        int score = (rand()%(int)CLR.size());
        for (CLR_t::const_iterator it=CLR.begin(); it!=CLR.end(); ++it)
        {
            if (count>=score)
            {
                Clrs.XPush (it->first);
                break;
            }
            count++;
        }
        if (Clrs.Size()==num) break;
        iter++;
    }
    if (Clrs.Size()!=num) throw new Fatal("Colors::GetRandom: failed with %d iterations",iter);
    size_t j = 0;
    for (size_t i=0; i<Num-num; ++i)
    {
        Clrs.Push (Clrs[j]);
        j++;
    }
}

char const * GetNext (size_t Idx)
{
    CLR_t::const_iterator p = CLR.begin();
    size_t k = Idx % CLR.size();
    for (size_t i=0; i<k; ++i) p++;
    return p->first.CStr();
}

char const * GetRandom ()
{
    size_t k = (rand()%(size_t)CLR.size());
    return GetNext (k);
}

int __init_colors__()
{
    // 191 colors
    Vec3_t c;
    c = 0.0, 0.0, 0.0;

    //  Whites
    c = 0.9804, 0.9216, 0.8431;  CLR["antique_white"]        = c;
    c = 0.9412, 1.0000, 1.0000;  CLR["azure"]                = c;
    c = 1.0000, 0.8941, 0.7686;  CLR["bisque"]               = c;
    c = 1.0000, 0.9216, 0.8039;  CLR["blanched_almond"]      = c;
    c = 1.0000, 0.9725, 0.8627;  CLR["cornsilk"]             = c;
    c = 0.9900, 0.9000, 0.7900;  CLR["eggshell"]             = c;
    c = 1.0000, 0.9804, 0.9412;  CLR["floral_white"]         = c;
    c = 0.8627, 0.8627, 0.8627;  CLR["gainsboro"]            = c;
    c = 0.9725, 0.9725, 1.0000;  CLR["ghost_white"]          = c;
    c = 0.9412, 1.0000, 0.9412;  CLR["honeydew"]             = c;
    c = 1.0000, 1.0000, 0.9412;  CLR["ivory"]                = c;
    c = 0.9020, 0.9020, 0.9804;  CLR["lavender"]             = c;
    c = 1.0000, 0.9412, 0.9608;  CLR["lavender_blush"]       = c;
    c = 1.0000, 0.9804, 0.8039;  CLR["lemon_chiffon"]        = c;
    c = 0.9804, 0.9412, 0.9020;  CLR["linen"]                = c;
    c = 0.9608, 1.0000, 0.9804;  CLR["mint_cream"]           = c;
    c = 1.0000, 0.8941, 0.8824;  CLR["misty_rose"]           = c;
    c = 1.0000, 0.8941, 0.7098;  CLR["moccasin"]             = c;
    c = 1.0000, 0.8706, 0.6784;  CLR["navajo_white"]         = c;
    c = 0.9922, 0.9608, 0.9020;  CLR["old_lace"]             = c;
    c = 1.0000, 0.9373, 0.8353;  CLR["papaya_whip"]          = c;
    c = 1.0000, 0.8549, 0.7255;  CLR["peach_puff"]           = c;
    c = 1.0000, 0.9608, 0.9333;  CLR["seashell"]             = c;
    c = 1.0000, 0.9804, 0.9804;  CLR["snow"]                 = c;
    c = 0.8471, 0.7490, 0.8471;  CLR["thistle"]              = c;
    c = 0.9900, 1.0000, 0.9400;  CLR["titanium_white"]       = c;
    c = 0.9608, 0.8706, 0.7020;  CLR["wheat"]                = c;
    c = 1.0000, 1.0000, 1.0000;  CLR["white"]                = c;
    c = 0.9608, 0.9608, 0.9608;  CLR["white_smoke"]          = c;
    c = 0.9900, 0.9700, 1.0000;  CLR["zinc_white"]           = c;

    // Greys
    c = 0.5000, 0.5400, 0.5300;  CLR["cold_grey"]            = c;
    c = 0.4118, 0.4118, 0.4118;  CLR["dim_grey"]             = c;
    c = 0.7529, 0.7529, 0.7529;  CLR["grey"]                 = c;
    c = 0.8275, 0.8275, 0.8275;  CLR["light_grey"]           = c;
    c = 0.4392, 0.5020, 0.5647;  CLR["slate_grey"]           = c;
    c = 0.1843, 0.3098, 0.3098;  CLR["slate_grey_dark"]      = c;
    c = 0.4667, 0.5333, 0.6000;  CLR["slate_grey_light"]     = c;
    c = 0.5000, 0.5000, 0.4100;  CLR["warm_grey"]            = c;

    // Blacks
    c = 0.0000, 0.0000, 0.0000;  CLR["black"]                = c;
    c = 0.1600, 0.1400, 0.1300;  CLR["ivory_black"]          = c;
    c = 0.1800, 0.2800, 0.2300;  CLR["lamp_black"]           = c;

    // Reds
    c = 0.8900, 0.1500, 0.2100;  CLR["alizarin_crimson"]     = c;
    c = 0.6100, 0.4000, 0.1200;  CLR["brick"]                = c;
    c = 0.8900, 0.0900, 0.0500;  CLR["cadmium_red_deep"]     = c;
    c = 1.0000, 0.4980, 0.3137;  CLR["coral"]                = c;
    c = 0.9412, 0.5020, 0.5020;  CLR["coral_light"]          = c;
    c = 1.0000, 0.0784, 0.5765;  CLR["deep_pink"]            = c;
    c = 0.8300, 0.2400, 0.1000;  CLR["english_red"]          = c;
    c = 0.6980, 0.1333, 0.1333;  CLR["firebrick"]            = c;
    c = 0.8900, 0.0700, 0.1900;  CLR["geranium_lake"]        = c;
    c = 1.0000, 0.4118, 0.7059;  CLR["hot_pink"]             = c;
    c = 0.6900, 0.0900, 0.1200;  CLR["indian_red"]           = c;
    c = 1.0000, 0.6275, 0.4784;  CLR["light_salmon"]         = c;
    c = 0.8900, 0.1800, 0.1900;  CLR["madder_lake_deep"]     = c;
    c = 0.6902, 0.1882, 0.3765;  CLR["maroon"]               = c;
    c = 1.0000, 0.7529, 0.7961;  CLR["pink"]                 = c;
    c = 1.0000, 0.7137, 0.7569;  CLR["pink_light"]           = c;
    c = 0.5300, 0.1500, 0.3400;  CLR["raspberry"]            = c;
    c = 1.0000, 0.0000, 0.0000;  CLR["red"]                  = c;
    c = 0.8900, 0.2100, 0.2200;  CLR["rose_madder"]          = c;
    c = 0.9804, 0.5020, 0.4471;  CLR["salmon"]               = c;
    c = 1.0000, 0.3882, 0.2784;  CLR["tomato"]               = c;
    c = 0.8300, 0.1000, 0.1200;  CLR["venetian_red"]         = c;

    // Browns
    c = 0.6400, 0.5800, 0.5000;  CLR["beige"]                = c;
    c = 0.5000, 0.1647, 0.1647;  CLR["brown"]                = c;
    c = 0.8600, 0.1600, 0.1600;  CLR["brown_madder"]         = c;
    c = 0.5300, 0.2600, 0.1200;  CLR["brown_ochre"]          = c;
    c = 0.8706, 0.7216, 0.5294;  CLR["burlywood"]            = c;
    c = 0.5400, 0.2100, 0.0600;  CLR["burnt_sienna"]         = c;
    c = 0.5400, 0.2000, 0.1400;  CLR["burnt_umber"]          = c;
    c = 0.8235, 0.4118, 0.1176;  CLR["chocolate"]            = c;
    c = 0.4500, 0.2400, 0.1000;  CLR["deep_ochre"]           = c;
    c = 1.0000, 0.4900, 0.2500;  CLR["flesh"]                = c;
    c = 1.0000, 0.3400, 0.1300;  CLR["flesh_ochre"]          = c;
    c = 0.7800, 0.4700, 0.1500;  CLR["gold_ochre"]           = c;
    c = 1.0000, 0.2400, 0.0500;  CLR["greenish_umber"]       = c;
    c = 0.9412, 0.9020, 0.5490;  CLR["khaki"]                = c;
    c = 0.7412, 0.7176, 0.4196;  CLR["khaki_dark"]           = c;
    c = 0.9608, 0.9608, 0.8627;  CLR["light_beige"]          = c;
    c = 0.8039, 0.5216, 0.2471;  CLR["peru"]                 = c;
    c = 0.7373, 0.5608, 0.5608;  CLR["rosy_brown"]           = c;
    c = 0.7800, 0.3800, 0.0800;  CLR["raw_sienna"]           = c;
    c = 0.4500, 0.2900, 0.0700;  CLR["raw_umber"]            = c;
    c = 0.3700, 0.1500, 0.0700;  CLR["sepia"]                = c;
    c = 0.6275, 0.3216, 0.1765;  CLR["sienna"]               = c;
    c = 0.5451, 0.2706, 0.0745;  CLR["saddle_brown"]         = c;
    c = 0.9569, 0.6431, 0.3765;  CLR["sandy_brown"]          = c;
    c = 0.8235, 0.7059, 0.5490;  CLR["tan"]                  = c;
    c = 0.3700, 0.1500, 0.0200;  CLR["van_dyke_brown"]       = c;

    // Oranges
    c = 1.0000, 0.3800, 0.0100;  CLR["cadmium_orange"]       = c;
    c = 1.0000, 0.0100, 0.0500;  CLR["cadmium_red_light"]    = c;
    c = 0.9300, 0.5700, 0.1300;  CLR["carrot"]               = c;
    c = 1.0000, 0.5490, 0.0000;  CLR["dark_orange"]          = c;
    c = 0.5900, 0.2700, 0.0800;  CLR["mars_orange"]          = c;
    c = 0.8900, 0.4400, 0.1000;  CLR["mars_yellow"]          = c;
    c = 1.0000, 0.5000, 0.0000;  CLR["orange"]               = c;
    c = 1.0000, 0.2706, 0.0000;  CLR["orange_red"]           = c;
    c = 0.8900, 0.5100, 0.0900;  CLR["yellow_ochre"]         = c;

    // Yellows
    c = 1.0000, 0.6600, 0.1400;  CLR["aureoline_yellow"]     = c;
    c = 0.8900, 0.8100, 0.3400;  CLR["banana"]               = c;
    c = 1.0000, 0.8900, 0.0100;  CLR["cadmium_lemon"]        = c;
    c = 1.0000, 0.6000, 0.0700;  CLR["cadmium_yellow"]       = c;
    c = 1.0000, 0.6900, 0.0600;  CLR["cadmium_yellow_light"] = c;
    c = 1.0000, 0.8431, 0.0000;  CLR["gold"]                 = c;
    c = 0.8549, 0.6471, 0.1255;  CLR["goldenrod"]            = c;
    c = 0.7216, 0.5255, 0.0431;  CLR["goldenrod_dark"]       = c;
    c = 0.9804, 0.9804, 0.8235;  CLR["goldenrod_light"]      = c;
    c = 0.9333, 0.9098, 0.6667;  CLR["goldenrod_pale"]       = c;
    c = 0.9333, 0.8667, 0.5098;  CLR["light_goldenrod"]      = c;
    c = 0.8900, 0.6600, 0.4100;  CLR["melon"]                = c;
    c = 1.0000, 0.6600, 0.0700;  CLR["naples_yellow_deep"]   = c;
    c = 1.0000, 1.0000, 0.0000;  CLR["yellow"]               = c;
    c = 1.0000, 1.0000, 0.8784;  CLR["yellow_light"]         = c;

    // Greens
    c = 0.4980, 1.0000, 0.0000;  CLR["chartreuse"]           = c;
    c = 0.4000, 0.5000, 0.0800;  CLR["chromeoxidegreen"]     = c;
    c = 0.3800, 0.7000, 0.1600;  CLR["cinnabar_green"]       = c;
    c = 0.2400, 0.5700, 0.2500;  CLR["cobalt_green"]         = c;
    c = 0.0000, 0.7900, 0.3400;  CLR["emerald_green"]        = c;
    c = 0.1333, 0.5451, 0.1333;  CLR["forest_green"]         = c;
    c = 0.0000, 1.0000, 0.0000;  CLR["green"]                = c;
    c = 0.0000, 0.3922, 0.0000;  CLR["green_dark"]           = c;
    c = 0.5961, 0.9843, 0.5961;  CLR["green_pale"]           = c;
    c = 0.6784, 1.0000, 0.1843;  CLR["green_yellow"]         = c;
    c = 0.4863, 0.9882, 0.0000;  CLR["lawn_green"]           = c;
    c = 0.1961, 0.8039, 0.1961;  CLR["lime_green"]           = c;
    c = 0.7400, 0.9900, 0.7900;  CLR["mint"]                 = c;
    c = 0.2300, 0.3700, 0.1700;  CLR["olive"]                = c;
    c = 0.4196, 0.5569, 0.1373;  CLR["olive_drab"]           = c;
    c = 0.3333, 0.4196, 0.1843;  CLR["olive_green_dark"]     = c;
    c = 0.0400, 0.7900, 0.1700;  CLR["permanent_green"]      = c;
    c = 0.1900, 0.5000, 0.0800;  CLR["sap_green"]            = c;
    c = 0.1804, 0.5451, 0.3412;  CLR["sea_green"]            = c;
    c = 0.5608, 0.7373, 0.5608;  CLR["sea_green_dark"]       = c;
    c = 0.2353, 0.7020, 0.4431;  CLR["sea_green_medium"]     = c;
    c = 0.1255, 0.6980, 0.6667;  CLR["sea_green_light"]      = c;
    c = 0.0000, 1.0000, 0.4980;  CLR["spring_green"]         = c;
    c = 0.0000, 0.9804, 0.6039;  CLR["spring_green_medium"]  = c;
    c = 0.2200, 0.3700, 0.0600;  CLR["terre_verte"]          = c;
    c = 0.4300, 1.0000, 0.4400;  CLR["viridian_light"]       = c;
    c = 0.6039, 0.8039, 0.1961;  CLR["yellow_green"]         = c;

    // Cyans
    c = 0.4980, 1.0000, 0.8314;  CLR["aquamarine"]           = c;
    c = 0.4000, 0.8039, 0.6667;  CLR["aquamarine_medium"]    = c;
    c = 0.0000, 1.0000, 1.0000;  CLR["cyan"]                 = c;
    c = 0.8784, 1.0000, 1.0000;  CLR["cyan_white"]           = c;
    c = 0.2510, 0.8784, 0.8157;  CLR["turquoise"]            = c;
    c = 0.0000, 0.8078, 0.8196;  CLR["turquoise_dark"]       = c;
    c = 0.2824, 0.8196, 0.8000;  CLR["turquoise_medium"]     = c;
    c = 0.6863, 0.9333, 0.9333;  CLR["turquoise_pale"]       = c;

    // Blues
    c = 0.9412, 0.9725, 1.0000;  CLR["alice_blue"]           = c;
    c = 0.0000, 0.0000, 1.0000;  CLR["blue"]                 = c;
    c = 0.6784, 0.8471, 0.9020;  CLR["blue_light"]           = c;
    c = 0.0000, 0.0000, 0.8039;  CLR["blue_medium"]          = c;
    c = 0.3725, 0.6196, 0.6275;  CLR["cadet"]                = c;
    c = 0.2400, 0.3500, 0.6700;  CLR["cobalt"]               = c;
    c = 0.3922, 0.5843, 0.9294;  CLR["cornflower"]           = c;
    c = 0.0200, 0.7200, 0.8000;  CLR["cerulean"]             = c;
    c = 0.1176, 0.5647, 1.0000;  CLR["dodger_blue"]          = c;
    c = 0.0300, 0.1800, 0.3300;  CLR["indigo"]               = c;
    c = 0.0100, 0.6600, 0.6200;  CLR["manganese_blue"]       = c;
    c = 0.0980, 0.0980, 0.4392;  CLR["midnight_blue"]        = c;
    c = 0.0000, 0.0000, 0.5020;  CLR["navy"]                 = c;
    c = 0.2000, 0.6300, 0.7900;  CLR["peacock"]              = c;
    c = 0.6902, 0.8784, 0.9020;  CLR["powder_blue"]          = c;
    c = 0.2549, 0.4118, 0.8824;  CLR["royal_blue"]           = c;
    c = 0.4157, 0.3529, 0.8039;  CLR["slate_blue"]           = c;
    c = 0.2824, 0.2392, 0.5451;  CLR["slate_blue_dark"]      = c;
    c = 0.5176, 0.4392, 1.0000;  CLR["slate_blue_light"]     = c;
    c = 0.4824, 0.4078, 0.9333;  CLR["slate_blue_medium"]    = c;
    c = 0.5294, 0.8078, 0.9216;  CLR["sky_blue"]             = c;
    c = 0.0000, 0.7490, 1.0000;  CLR["sky_blue_deep"]        = c;
    c = 0.5294, 0.8078, 0.9804;  CLR["sky_blue_light"]       = c;
    c = 0.2745, 0.5098, 0.7059;  CLR["steel_blue"]           = c;
    c = 0.6902, 0.7686, 0.8706;  CLR["steel_blue_light"]     = c;
    c = 0.0000, 0.7800, 0.5500;  CLR["turquoise_blue"]       = c;
    c = 0.0700, 0.0400, 0.5600;  CLR["ultramarine"]          = c;

    // Magenta
    c = 0.5412, 0.1686, 0.8863;  CLR["blue_violet"]          = c;
    c = 0.5700, 0.1300, 0.6200;  CLR["cobalt_violet_deep"]   = c;
    c = 1.0000, 0.0000, 1.0000;  CLR["magenta"]              = c;
    c = 0.8549, 0.4392, 0.8392;  CLR["orchid"]               = c;
    c = 0.6000, 0.1961, 0.8000;  CLR["orchid_dark"]          = c;
    c = 0.7294, 0.3333, 0.8275;  CLR["orchid_medium"]        = c;
    c = 0.8600, 0.1500, 0.2700;  CLR["permanent_red_violet"] = c;
    c = 0.8667, 0.6275, 0.8667;  CLR["plum"]                 = c;
    c = 0.6275, 0.1255, 0.9412;  CLR["purple"]               = c;
    c = 0.5765, 0.4392, 0.8588;  CLR["purple_medium"]        = c;
    c = 0.3600, 0.1400, 0.4300;  CLR["ultramarine_violet"]   = c;
    c = 0.5600, 0.3700, 0.6000;  CLR["violet"]               = c;
    c = 0.5804, 0.0000, 0.8275;  CLR["violet_dark"]          = c;
    c = 0.8157, 0.1255, 0.5647;  CLR["violet_red"]           = c;
    c = 0.7804, 0.0824, 0.5216;  CLR["violet_red_medium"]    = c;
    c = 0.8588, 0.4392, 0.5765;  CLR["violet_red_pale"]      = c;

    return 0;
}

int __colors_dummy__ = __init_colors__();

}; // namespace Colors

#endif // MECHSYS_COLORS_H
