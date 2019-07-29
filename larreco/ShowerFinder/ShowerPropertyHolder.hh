//###################################################################
//### Name:        ShowerPropertyHolder                           ### 
//### Author:      Dominic Barker                                 ###
//### Date:        15.07.19                                       ###
//### Description: Class to holder the standard shower property   ###
//###              information. Used in SBNShower_modle and       ###
//###              corresponding tools.                           ###
//###################################################################

#ifndef ShowerPropertyHolder_HH
#define ShowerPropertyHolder_HH

//C++ Inlcudes 
#include <iostream> 
#include <map>
#include <string> 
#include <memory>

namespace reco {
  namespace shower {
    template <class T> class ShowerProperty;
    class ShowerPropertyBase;
    class ShowerPropertyHolder;
  }
}

class reco::shower::ShowerPropertyBase {


public:

  virtual ~ShowerPropertyBase() noexcept = default;
  
  bool CheckShowerProperty(){
    if(propertyPtr) return true;
    else return false;
  }

  bool CheckShowerPropertyError(){
    if(propertyErrPtr) return true;
    else return false;
  }

  void Clear(){
    propertyPtr    = 0;
    propertyErrPtr = 0;
  }


protected:

  int propertyPtr;
  int propertyErrPtr; 
  
};

//This is a template class which holds a shower property. This holds any object e.g. std::vector<double>, double, TVector3 
//and holds various information which helps the showerproperty holder access the elements. A user should not require any part 
//of this class.  
template <class T>
class reco::shower::ShowerProperty : public reco::shower::ShowerPropertyBase {

public: 
  
  ShowerProperty(){
    propertyPtr    = 0;
    propertyErrPtr = 0;
  }

  void SetShowerProperty(T& Property){
    property = Property;
    propertyPtr = 1;
  }

  void SetShowerPropertryError(T& PropertyError){
    propertyErr = PropertyError;
    propertyErrPtr = 1;
  }
  
  int GetShowerProperty(T& Property){
    if(propertyPtr){
      Property = property;
      return 0;
    }
    else{
      return 1;
    }
  }

  int GetShowerPropertyError(T& PropertyErr){
    if(propertyErrPtr){
      PropertyErr = propertyErr;
      return 0;
    }
    else{
      return 1;
    }
  }
    
private:
  T   property; 
  T   propertyErr;
};


//Class to holder all the reco::shower::ShowerProperty objects. This is essentially a map from a string the object so people can 
//add an object in a tool and get it back later. 
class reco::shower::ShowerPropertyHolder{

public:
  
  //Getter function for accessing the shower property e..g the direction ShowerPropertyHolder.GetProperty("MyShowerValue"); The name is used access the value and precise names are required for a complete shower in sbnshower: ShowerStartPosition, ShowerDirection, ShowerEnergy ,ShowerdEdx.
  template <class T>
  int GetProperty(std::string Name, T& Property){
    if(showerproperties.find(Name) == showerproperties.end()){return 1;}
    //    reco::shower::ShowerPropertyBase * showerprop_base = showerproperties[Name];
    reco::shower::ShowerProperty<T> *showerprop = dynamic_cast<reco::shower::ShowerProperty<T> *>(showerproperties[Name].get());
    showerprop->GetShowerProperty(Property);
    return 0;
  }

  //Getter function for accessing the shower property error e.g the direction ShowerPropertyHolder.GetProperty("MyShowerValue");
  template <class T>
  int GetPropertyError(std::string Name, T& PropertyErr){
    if(showerproperties.find(Name) == showerproperties.end()){return 1;}
    reco::shower::ShowerProperty<T> *showerprop = dynamic_cast<reco::shower::ShowerProperty<T> *>(showerproperties[Name].get());
    return showerprop->GetShowerPropertyError(PropertyErr);
  }

  //This sets the value of the property. Just give a name and a object
  //e.g. TVector3 ShowerPropertyHolder.SetProperty((TVector3) StartPosition, "StartPosition");
  template <class T>
  void SetProperty(T& propertyval, std::string Name){

    if(showerproperties.find(Name) != showerproperties.end()){
      reco::shower::ShowerProperty<T>* showerprop = dynamic_cast<reco::shower::ShowerProperty<T> *>(showerproperties[Name].get());
      showerprop->SetShowerProperty(propertyval);
    return;
    }
    else{
      reco::shower::ShowerProperty<T> * ptr = new reco::shower::ShowerProperty<T>();
      showerproperties[Name] =  std::unique_ptr<reco::shower::ShowerPropertyBase>(ptr);
      reco::shower::ShowerProperty<T>* showerprop = dynamic_cast<reco::shower::ShowerProperty<T> *>(showerproperties[Name].get());
      showerprop->SetShowerProperty(propertyval);
      return;
    }
  }

  //This sets the value of the property error. Just give a name and a object
  //e.g. TVector3 ShowerPropertyHolder.SetPropertyError((TVector3) StartPosition, "StartPoision");
  template <class T>
  void SetPropertyError(T propertyval, std::string Name){

    if(showerproperties.find(Name) == showerproperties.end()){
      ((reco::shower::ShowerProperty<T>) *showerproperties[Name]).SetShowerPropertyError(propertyval);
      return;
    }
    else{
      showerproperties[Name] = std::make_unique<reco::shower::ShowerProperty<T> >();
      ((reco::shower::ShowerProperty<T>) *showerproperties[Name]).SetShowerPropertyError(propertyval);
      return;
    }
  }

  //Check that a property is filled 
  bool CheckProperty(std::string Name){
    if(showerproperties.find(Name) == showerproperties.end()){return false;}
    return showerproperties[Name]->CheckShowerProperty();
  }
  
  //Check that the property error is filled 
  bool CheckPropertyError(std::string Name){
    if(showerproperties.find(Name) == showerproperties.end()){return false;}
    return showerproperties[Name]->CheckShowerPropertyError();
  }

  //Check All the properties
  bool CheckAllProperties(){
    bool checked = true;
    for(auto const& showerprop: showerproperties){
     checked *= showerprop.second->CheckShowerProperty();
    }
    return checked;
  }
  
  //Check All the properties errors
  bool CheckAllPropertiesErrors(){
    bool checked = true;
    for(auto const& showerprop: showerproperties){
      checked *= (showerprop.second)->CheckShowerPropertyError();
    }
    return checked;
  }

  
  //Clear Fucntion
  void ClearProperty(std::string Name){
    if(showerproperties.find(Name) == showerproperties.end()){return;}
    return showerproperties[Name]->Clear();
  }

  //Clear all the shower propertys
  void ClearAll(){
    for(auto const& showerprop: showerproperties){
      (showerprop.second)->Clear();
    }
  }

private:

  //Storage for all the shower properties.
  std::map<std::string,std::unique_ptr<reco::shower::ShowerPropertyBase> > showerproperties;

};


#endif
